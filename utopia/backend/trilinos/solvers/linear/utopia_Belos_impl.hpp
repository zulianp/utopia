#ifndef UTOPIA_BELOS_IMPL_HPP
#define UTOPIA_BELOS_IMPL_HPP

#include "utopia_Belos_solver.hpp"

#include "utopia_make_unique.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp> //TODO remove from here


#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>

//TODO remove from here
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>


#ifdef WITH_TRILINOS_MUELU
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#else
#warning "Trilinos was not built with MueLu support. AMG cannot be used as a preconditioner for the Belos solver."
#endif //WITH_TRILINOS_MUELU


#ifdef WITH_TRILINOS_IFPACK2
#include <Ifpack2_Factory.hpp>
#else
#warning "Trilinos was not built with Ifpack2 support. Direct preconditioners cannot be used with the Belos solver."
#endif //WITH_TRILINOS_IFPACK2


namespace utopia {

    template <typename Matrix, typename Vector>
    class BelosSolver<Matrix, Vector, TRILINOS>::Impl {
    public:
        typedef double ST;

        typedef Tpetra::Operator<>::scalar_type SC;
        typedef Tpetra::Operator<SC>::local_ordinal_type LO;
        typedef Tpetra::Operator<SC, LO>::global_ordinal_type GO;

        typedef Kokkos::Compat::KokkosSerialWrapperNode serial_node;

#ifdef  KOKKOS_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode cuda_node;
        typedef cuda_node NT;
#elif defined   KOKKOS_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode openmp_node;
        typedef Kokkos::Compat::KokkosThreadsWrapperNode thread_node;
        typedef openmp_node NT;
#else
        typedef serial_node NT;
#endif

        typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
        typedef Tpetra::Operator<SC, LO, GO, NT> OP;

        typedef Belos::LinearProblem<SC, MV, OP> problem_type;
        typedef Belos::SolverManager<SC, MV, OP> solver_type;

#ifdef WITH_TRILINOS_IFPACK2
        typedef Ifpack2::Preconditioner<SC, LO, GO, NT> ifpack_prec_type;
#endif //WITH_TRILINOS_IFPACK2

#ifdef WITH_TRILINOS_MUELU
        typedef MueLu::TpetraOperator<SC, LO, GO, NT> muelu_prec_type;
#endif

        typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;
        typedef Tpetra::CrsMatrix<SC, LO, GO, NT> matrix_type;

        Teuchos::RCP<problem_type> linear_problem;
        Teuchos::RCP<Teuchos::ParameterList> param_list;
        //  auto& utopiaPL;// impl_->param_list->sublist("UTOPIA", true);
        Teuchos::RCP<solver_type> belos_solver;
        Belos::SolverFactory<SC, MV, OP> belos_factory;

        //preconditioner
#ifdef WITH_TRILINOS_IFPACK2
        Teuchos::RCP<ifpack_prec_type> M_ifpack;
#endif //WITH_TRILINOS_IFPACK2

#ifdef WITH_TRILINOS_MUELU
        Teuchos::RCP<muelu_prec_type> M_muelu;
#endif //WITH_TRILINOS_MUELU

    };

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver(Parameters params)
    : impl_(make_unique<Impl>())
    {
        set_parameters(params);
    }

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver(const BelosSolver &other)
    : impl_(make_unique<Impl>(*other.impl_)) {
        //FIXME
    }

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::~BelosSolver() {}

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver() : impl_(make_unique<Impl>()) {}

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op,
                                                       const std::shared_ptr<const Matrix> &prec)
    {
        PreconditionedSolver::update(op, prec);
        // set_problem(*op);
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op)
    {
        PreconditionedSolver::update(op);
        // set_problem(*op);
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::apply(const Vector &rhs, Vector &lhs) {

        impl_->linear_problem = Teuchos::rcp(
                                             new typename Impl::problem_type(
                                                                             this->get_operator()->implementation().implementation_ptr(),
                                                                             lhs.implementation().implementation_ptr(),
                                                                             rhs.implementation().implementation_ptr()
                                                                             )
                                             );

        set_problem();

        assert(!(impl_->belos_solver.is_null()));
        impl_->belos_solver->solve();
        return true;
    }

    template <typename Matrix, typename Vector>
    int BelosSolver<Matrix, Vector, TRILINOS>::get_num_iter() const { return impl_->belos_solver->getNumIters(); }


    template <typename Matrix, typename Vector>
    double BelosSolver<Matrix, Vector, TRILINOS>::achieved_tol() const { return impl_->belos_solver->achievedTol(); }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
    {
        set_preconditioner(); //(A) //FIXME
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const Matrix &precond)
    {
        bool direct_solver = impl_->param_list->sublist("UTOPIA", true).template get<bool>("Direct Preconditioner", false);
        std::string dir_prec_type = impl_->param_list->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");
        if ( direct_solver ) {
#ifdef WITH_TRILINOS_IFPACK2
            impl_->M_ifpack = Ifpack2::Factory::create<typename Impl::matrix_type>(dir_prec_type, precond.implementation().implementation_ptr());
            assert(!impl_->M_ifpack.is_null());
            impl_->M_ifpack->setParameters(impl_->param_list->sublist(dir_prec_type, false));
            impl_->M_ifpack->initialize();
            impl_->M_ifpack->compute();
            impl_->linear_problem->setRightPrec(impl_->M_ifpack);
#else  // WITH_TRILINOS_IFPACK2
          std::cerr << "Cannot use a Direct Preconditioner with the BelosSolver, since Trilinos was not built with Ifpack2 support!" << std::endl;
#endif // WITH_TRILINOS_IFPACK2
        } else {
#ifdef WITH_TRILINOS_MUELU
            // Multigrid Hierarchy
            impl_->M_muelu = MueLu::CreateTpetraPreconditioner((
                                                                Teuchos::RCP<typename Impl::OP>) precond.implementation().implementation_ptr(),
                                                               impl_->param_list->sublist("MueLu", false)
                                                               );

            assert(!impl_->M_muelu.is_null());
            impl_->linear_problem->setRightPrec(impl_->M_muelu);
#else
            std::cerr << "Cannot use MueLu as preconditioner since Trilinos was not built with MueLu support." << std::endl;
#endif //WITH_TRILINOS_MUELU
        }
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_parameters(const Parameters params)
    {
        if(!params.param_file_name().empty()) {
            try {
                impl_->param_list = Teuchos::getParametersFromXmlFile(params.param_file_name());
            } catch(const std::exception &ex) {
                std::cerr << ex.what() << std::endl;
                assert(false);
                abort();
            }
        } else {
            //use default paramlist
        }
    }

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS> * BelosSolver<Matrix, Vector, TRILINOS>::clone() const
    {
        return new BelosSolver(*this);
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::smooth(const Vector &rhs, Vector &x)
    {
        return false;
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::set_problem()
    {
        impl_->linear_problem->setProblem();
        auto sol_type = impl_->param_list->get("Solver Type", "CG");
        auto belos_params = Teuchos::sublist(impl_->param_list, sol_type, false);
        impl_->belos_solver = impl_->belos_factory.create( sol_type, belos_params); //to change it to have the specialization
        impl_->belos_solver->setProblem(impl_->linear_problem);
        if (this->verbose()) { impl_->belos_solver->getCurrentParameters()->print(); }
        set_preconditioner();
        return true;
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::set_problem(Matrix &A)
    {
        impl_->linear_problem->setProblem();
        auto sol_type = impl_->param_list->get("Solver Type", "CG");
        auto belos_params = Teuchos::sublist(impl_->param_list, sol_type, false);
        impl_->belos_solver = impl_->belos_factory.create( sol_type, belos_params); //to change it to have the specialization
        set_preconditioner(); //(A);
        impl_->belos_solver->setProblem(impl_->linear_problem);
        if (this->verbose()) { impl_->belos_solver->getCurrentParameters()->print(); }
        return true;
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner()//const std::shared_ptr<Preconditioner> &precond)
    {
        bool direct_solver = impl_->param_list->sublist("UTOPIA", true).template get<bool>("Direct Preconditioner", false);
        std::string dir_prec_type = impl_->param_list->sublist("UTOPIA", true).get("Ifpack2 Preconditioner", "prec_type_unset");

        if ( direct_solver )
        {
#ifdef WITH_TRILINOS_IFPACK2
            impl_->M_ifpack = Ifpack2::Factory::create<typename Impl::matrix_type>(dir_prec_type, this->get_operator()->implementation().implementation_ptr());
            assert(!impl_->M_ifpack.is_null());
            impl_->M_ifpack->setParameters(impl_->param_list->sublist(dir_prec_type, false));
            impl_->M_ifpack->initialize();
            impl_->M_ifpack->compute();
            impl_->linear_problem->setRightPrec(impl_->M_ifpack);
#else  //WITH_TRILINOS_IFPACK2
          std::cerr << "Cannot use a Direct Preconditioner with the BelosSolver, since Trilinos was not built with Ifpack2 support!" << std::endl;
#endif //WITH_TRILINOS_IFPACK2
        } else {
#ifdef WITH_TRILINOS_MUELU
            // Multigrid Hierarchy
            //                M_muelu = MueLu::CreateTpetraPreconditioner((Teuchos::RCP<OP>)precond->implementation().implementation_ptr(),    //TODO
            //                                                            impl_->param_list->sublist("MueLu", false));
            assert(!impl_->M_muelu.is_null());
            impl_->linear_problem->setRightPrec(impl_->M_muelu);
#else  // WITH_TRILINOS_MUELU
            std::cerr << "Cannot use MueLu as preconditioner since Trilinos was not built with MueLu support." << std::endl;
#endif //WITH_TRILINOS_MUELU
        }
    }

}  // namespace utopia

#endif //UTOPIA_BELOS_IMPL_HPP
