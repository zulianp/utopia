#ifndef UTOPIA_BELOS_IMPL_HPP
#define UTOPIA_BELOS_IMPL_HPP

#include "utopia_Belos_solver.hpp"

#include "utopia_Wrapper.hpp"
#include "utopia_make_unique.hpp"

#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

// TODO remove from here
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>

#ifdef WITH_TRILINOS_MUELU
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#else
#warning "Trilinos was not built with MueLu support. AMG cannot be used as a preconditioner for the Belos solver."
#endif  // WITH_TRILINOS_MUELU

#ifdef WITH_TRILINOS_IFPACK2
#include <Ifpack2_Factory.hpp>
#else
#warning "Trilinos was not built with Ifpack2 support. Direct preconditioners cannot be used with the Belos solver."
#endif  // WITH_TRILINOS_IFPACK2

namespace utopia {

    template <typename Matrix, typename Vector>
    class BelosSolver<Matrix, Vector, TRILINOS>::Impl {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using LocalSizeType = typename Traits<Vector>::LocalSizeType;
        using Node = typename Traits<Vector>::Node;

        using VectorType = typename Vector::VectorType;
        using MultiVectorType = typename Vector::MultiVectorType;
        using CrsMatrixType = typename Matrix::CrsMatrixType;
        using OperatorType = Tpetra::Operator<Scalar, LocalSizeType, SizeType, Node>;

        using ProblemType = Belos::LinearProblem<Scalar, MultiVectorType, OperatorType>;
        using SolverType = Belos::SolverManager<Scalar, MultiVectorType, OperatorType>;

#ifdef WITH_TRILINOS_IFPACK2
        using IfPack2PrecType = Ifpack2::Preconditioner<Scalar, LocalSizeType, SizeType, Node>;
#endif  // WITH_TRILINOS_IFPACK2

#ifdef WITH_TRILINOS_MUELU
        using MueLuPrecType = MueLu::TpetraOperator<Scalar, LocalSizeType, SizeType, Node>;
#endif

        Teuchos::RCP<ProblemType> linear_problem;
        Teuchos::RCP<Teuchos::ParameterList> param_list;
        Teuchos::RCP<SolverType> belos_solver;
        Belos::SolverFactory<Scalar, MultiVectorType, OperatorType> belos_factory;

        // preconditioner
#ifdef WITH_TRILINOS_IFPACK2
        Teuchos::RCP<IfPack2PrecType> ifpack2_prec_;
#endif  // WITH_TRILINOS_IFPACK2

#ifdef WITH_TRILINOS_MUELU
        Teuchos::RCP<MueLuPrecType> muelu_prec_;
#endif  // WITH_TRILINOS_MUELU
    };

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver(const BelosSolver &other)
        : PreconditionedSolver(other), impl_(utopia::make_unique<Impl>(*other.impl_)) {
        // FIXME
    }

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::~BelosSolver() {}

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver() : impl_(utopia::make_unique<Impl>()) {}

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op,
                                                       const std::shared_ptr<const Matrix> &prec) {
        PreconditionedSolver::update(op, prec);
        // set_problem(*op);
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op) {
        PreconditionedSolver::update(op);
        // set_problem(*op);
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::apply(const Vector &rhs, Vector &lhs) {
        impl_->linear_problem =
            Teuchos::rcp(new typename Impl::ProblemType(raw_type(*this->get_operator()), raw_type(lhs), raw_type(rhs)));

        set_problem();

        assert(!(impl_->belos_solver.is_null()));
        impl_->belos_solver->solve();
        return true;
    }

    template <typename Matrix, typename Vector>
    int BelosSolver<Matrix, Vector, TRILINOS>::get_num_iter() const {
        return impl_->belos_solver->getNumIters();
    }

    template <typename Matrix, typename Vector>
    double BelosSolver<Matrix, Vector, TRILINOS>::achieved_tol() const {
        return impl_->belos_solver->achievedTol();
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const std::shared_ptr<Preconditioner> &) {
        assert(false && "IMPLEMENT ME");
        set_preconditioner();  //(A) //FIXME
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const Matrix &precond) {
        bool direct_solver =
            impl_->param_list->sublist("UTOPIA", true).template get<bool>("Direct Preconditioner", false);
        std::string dir_prec_type =
            impl_->param_list->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
        if (direct_solver) {
#ifdef WITH_TRILINOS_IFPACK2
            impl_->ifpack2_prec_ =
                Ifpack2::Factory::create<typename Impl::CrsMatrixType>(dir_prec_type, raw_type(precond));
            assert(!impl_->ifpack2_prec_.is_null());
            impl_->ifpack2_prec_->setParameters(impl_->param_list->sublist(dir_prec_type, false));
            impl_->ifpack2_prec_->initialize();
            impl_->ifpack2_prec_->compute();
            std::string preconditioner_type =
                impl_->param_list->sublist("UTOPIA", true).get("Preconditioner Type", "right");
            std::transform(preconditioner_type.begin(),
                           preconditioner_type.end(),
                           preconditioner_type.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (preconditioner_type == "left") {
                impl_->linear_problem->setLeftPrec(impl_->ifpack2_prec_);
            } else {
                impl_->linear_problem->setRightPrec(impl_->ifpack2_prec_);
            }
#else   // WITH_TRILINOS_IFPACK2
            std::cerr << "Cannot use a Direct Preconditioner with the BelosSolver, since Trilinos was not built with "
                         "Ifpack2 support!"
                      << std::endl;
#endif  // WITH_TRILINOS_IFPACK2
        } else {
#ifdef WITH_TRILINOS_MUELU
            // Multigrid Hierarchy
            impl_->muelu_prec_ =
                MueLu::CreateTpetraPreconditioner((Teuchos::RCP<typename Impl::OperatorType>)raw_type(precond),
                                                  impl_->param_list->sublist("MueLu", false));

            assert(!impl_->muelu_prec_.is_null());
            std::string preconditioner_type =
                impl_->param_list->sublist("UTOPIA", true).get("Preconditioner Type", "right");
            std::transform(preconditioner_type.begin(),
                           preconditioner_type.end(),
                           preconditioner_type.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (preconditioner_type == "left") {
                impl_->linear_problem->setLeftPrec(impl_->muelu_prec_);
            } else {
                impl_->linear_problem->setRightPrec(impl_->muelu_prec_);
            }
#else
            std::cerr << "Cannot use MueLu as preconditioner since Trilinos was not built with MueLu support."
                      << std::endl;
#endif  // WITH_TRILINOS_MUELU
        }
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::read_xml(const std::string &path) {
        if (!path.empty()) {
            try {
                impl_->param_list = Teuchos::getParametersFromXmlFile(path);
            } catch (const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
                abort();
            }
        } else {
            // use default paramlist
        }
    }

    // read an utopia file and convert in a trilinos list
    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::read(Input &in) {
        PreconditionedSolver::read(in);

        // TODO
        std::string exotic = "";
        in.get("exotic", exotic);

        // if (!exotic.empty()) {
        // }

        if (impl_->param_list.is_null()) {
            impl_->param_list = Teuchos::parameterList();
        }

        impl_->param_list->set("Relative tolerance", this->rtol(), "CG");
        impl_->param_list->set("S tolerance", this->stol(), "CG");
        impl_->param_list->set("A tolerance", this->atol(), "CG");
        impl_->param_list->set("Maximum iteration", this->max_it(), "CG");
        impl_->param_list->set("Verbose", this->verbose(), "CG");
        // auto in = open_istream(const Path &path);
    }

    // available parameters
    // TODO print setted parameters??
    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::print_usage(std::ostream &os) const {
        PreconditionedSolver::print_usage(os);
        // TODO
        // m_utopia_warning_once("not implemented");
    }

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS> *BelosSolver<Matrix, Vector, TRILINOS>::clone() const {
        return new BelosSolver(*this);
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::set_problem() {
        impl_->linear_problem->setProblem();
        auto sol_type = impl_->param_list->get("Solver Type", "CG");
        auto belos_params = Teuchos::sublist(impl_->param_list, sol_type, false);
        impl_->belos_solver =
            impl_->belos_factory.create(sol_type, belos_params);  // to change it to have the specialization
        impl_->belos_solver->setProblem(impl_->linear_problem);
        if (this->verbose()) {
            impl_->belos_solver->getCurrentParameters()->print();
        }
        set_preconditioner();
        return true;
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::set_problem(Matrix &) {
        std::cerr << "[Warning] matrix parameter ignored" << std::endl;

        impl_->linear_problem->setProblem();
        auto sol_type = impl_->param_list->get("Solver Type", "CG");
        auto belos_params = Teuchos::sublist(impl_->param_list, sol_type, false);
        impl_->belos_solver =
            impl_->belos_factory.create(sol_type, belos_params);  // to change it to have the specialization
        set_preconditioner();                                     //(A);
        impl_->belos_solver->setProblem(impl_->linear_problem);
        if (this->verbose()) {
            impl_->belos_solver->getCurrentParameters()->print();
        }
        return true;
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner()  // const std::shared_ptr<Preconditioner> &precond)
    {
        bool direct_solver =
            impl_->param_list->sublist("UTOPIA", true).template get<bool>("Direct Preconditioner", false);
        std::string dir_prec_type =
            impl_->param_list->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");

        if (direct_solver) {
#ifdef WITH_TRILINOS_IFPACK2
            impl_->ifpack2_prec_ =
                Ifpack2::Factory::create<typename Impl::CrsMatrixType>(dir_prec_type, raw_type(*this->get_operator()));
            assert(!impl_->ifpack2_prec_.is_null());
            impl_->ifpack2_prec_->setParameters(impl_->param_list->sublist(dir_prec_type, false));
            impl_->ifpack2_prec_->initialize();
            impl_->ifpack2_prec_->compute();
            std::string preconditioner_type =
                impl_->param_list->sublist("UTOPIA", true).get("Preconditioner Type", "right");
            // TODO to move to input validation phase
            std::transform(preconditioner_type.begin(),
                           preconditioner_type.end(),
                           preconditioner_type.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (preconditioner_type == "left") {
                impl_->linear_problem->setLeftPrec(impl_->ifpack2_prec_);
            } else {
                impl_->linear_problem->setRightPrec(impl_->ifpack2_prec_);
            }
#else   // WITH_TRILINOS_IFPACK2
            std::cerr << "Cannot use a Direct Preconditioner with the BelosSolver, since Trilinos was not built with "
                         "Ifpack2 support!"
                      << std::endl;
#endif  // WITH_TRILINOS_IFPACK2
        } else {
#ifdef WITH_TRILINOS_MUELU
            // Multigrid Hierarchy
            impl_->muelu_prec_ = MueLu::CreateTpetraPreconditioner(
                (Teuchos::RCP<typename Impl::OperatorType>)raw_type(*this->get_operator()),
                impl_->param_list->sublist("MueLu", false));
            assert(!impl_->muelu_prec_.is_null());
            std::string preconditioner_type =
                impl_->param_list->sublist("UTOPIA", true).get("Preconditioner Type", "right");
            // TODO to move to input validation phase
            std::transform(preconditioner_type.begin(),
                           preconditioner_type.end(),
                           preconditioner_type.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (preconditioner_type == "left") {
                impl_->linear_problem->setLeftPrec(impl_->muelu_prec_);
            } else {
                impl_->linear_problem->setRightPrec(impl_->muelu_prec_);
            }
#else   // WITH_TRILINOS_MUELU
            std::cerr << "Cannot use MueLu as preconditioner since Trilinos was not built with MueLu support."
                      << std::endl;
#endif  // WITH_TRILINOS_MUELU
        }
    }

}  // namespace utopia

#endif  // UTOPIA_BELOS_IMPL_HPP
