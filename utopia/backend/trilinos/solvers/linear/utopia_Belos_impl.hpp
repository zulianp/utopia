#ifndef UTOPIA_BELOS_IMPL_HPP
#define UTOPIA_BELOS_IMPL_HPP

#include "utopia_Base.hpp"
#include "utopia_Belos_solver.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_make_unique.hpp"

#ifdef UTOPIA_ENABLE_TRILINOS_BELOS
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

#ifdef UTOPIA_ENABLE_TRILINOS_MUELU
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#else
#warning "Trilinos was built without MueLu support. AMG cannot be used as a preconditioner for the Belos solver."
#endif  // UTOPIA_ENABLE_TRILINOS_MUELU

#ifdef UTOPIA_ENABLE_TRILINOS_IFPACK2
#include <Ifpack2_Factory.hpp>
#else
#warning "Trilinos was built without Ifpack2 support. Direct preconditioners cannot be used with the Belos solver."
#endif  // UTOPIA_ENABLE_TRILINOS_IFPACK2

#include "utopia_Tpetra_Operator.hpp"

namespace utopia {

    template <typename Matrix, typename Vector>
    class BelosSolver<Matrix, Vector, TRILINOS>::Impl {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using LocalSizeType = typename Traits<Vector>::LocalSizeType;
        using Node = typename Traits<Vector>::Node;

        using MultiVectorType = typename Vector::MultiVectorType;
        using CrsMatrixType = typename Matrix::CrsMatrixType;
        using OperatorType = Tpetra::Operator<Scalar, LocalSizeType, SizeType, Node>;

        using ProblemType = Belos::LinearProblem<Scalar, MultiVectorType, OperatorType>;
        using SolverType = Belos::SolverManager<Scalar, MultiVectorType, OperatorType>;

#ifdef UTOPIA_ENABLE_TRILINOS_IFPACK2
        using IfPack2PrecType = Ifpack2::Preconditioner<Scalar, LocalSizeType, SizeType, Node>;
#endif  // UTOPIA_ENABLE_TRILINOS_IFPACK2

#ifdef UTOPIA_ENABLE_TRILINOS_MUELU
        using MueLuPrecType = MueLu::TpetraOperator<Scalar, LocalSizeType, SizeType, Node>;
#endif  // UTOPIA_ENABLE_TRILINOS_MUELU

        Teuchos::RCP<ProblemType> linear_problem;
        Teuchos::RCP<Teuchos::ParameterList> param_list;
        Teuchos::RCP<SolverType> belos_solver;
        Belos::SolverFactory<Scalar, MultiVectorType, OperatorType> belos_factory;

        Impl() : param_list(Teuchos::parameterList()) {}
    };

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver(const BelosSolver &other)
        : PreconditionedSolverInterface<Vector>(other),
          Super(other),
          solver_type_(other.solver_type_),
          solver_params_(other.solver_params_),
          impl_(utopia::make_unique<Impl>(*other.impl_)) {}

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::BelosSolver(const std::string &type, std::initializer_list<ParamKey> params)
        : solver_type_(type), solver_params_(params), impl_(utopia::make_unique<Impl>()) {}

    template <typename Matrix, typename Vector>
    BelosSolver<Matrix, Vector, TRILINOS>::~BelosSolver() {}

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op,
                                                       const std::shared_ptr<const Matrix> &prec) {
        Super::update(op, prec);
        // set_problem(*op);
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::update(const std::shared_ptr<const Matrix> &op) {
        Super::update(op);
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::apply(const Vector &rhs, Vector &lhs) {
        UTOPIA_TRACE_SCOPE("BelosSolver::apply(...)");
        if (empty(lhs)) {
            lhs.zeros(layout(rhs));
        }

        assert(!rhs.has_nan_or_inf());

        impl_->linear_problem =
            Teuchos::rcp(new typename Impl::ProblemType(raw_type(*this->get_operator()), raw_type(lhs), raw_type(rhs)));

        set_programmatic_config();
        set_problem();
        set_preconditioner(this->get_operator());

        assert(!(impl_->belos_solver.is_null()));
        impl_->belos_solver->solve();
        return true;
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::print_usage(std::ostream &os) const {
        Super::print_usage(os);
        for (const auto &param : solver_params_) {
            this->print_param_usage(
                os, params_[param].name, params_[param].type, params_[param].description, params_[param].default_value);
        }
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::read(Input &in) {
        Super::read(in);
        set_programmatic_config();

        // The solver_params sublist is created by the call to set_programmatic_config
        constexpr bool must_exist{true};
        auto solver_params = Teuchos::sublist(impl_->param_list, solver_type_, must_exist);

        // Set Belos-specific parameters
        if (this->verbose()) {
            in.get(params_[ParamKey::OUTPUT_FREQ].name, output_frequency_);
            solver_params->set(params_[ParamKey::OUTPUT_FREQ].belos_name, output_frequency_);
        }
        if (solver_params_.find(ParamKey::BLOCK_SIZE) != solver_params_.end()) {
            in.get(params_[ParamKey::BLOCK_SIZE].name, block_size_);
            solver_params->set(params_[ParamKey::BLOCK_SIZE].belos_name, block_size_);
        }
        if (solver_params_.find(ParamKey::MAX_RESTARTS) != solver_params_.end()) {
            in.get(params_[ParamKey::MAX_RESTARTS].name, max_restarts_);
            solver_params->set(params_[ParamKey::MAX_RESTARTS].belos_name, max_restarts_);
        }
        if (solver_params_.find(ParamKey::NUM_BLOCKS) != solver_params_.end()) {
            in.get(params_[ParamKey::NUM_BLOCKS].name, num_blocks_);
            solver_params->set(params_[ParamKey::NUM_BLOCKS].belos_name, num_blocks_);
        }
        if (solver_params_.find(ParamKey::ORTHOGONALIZATION) != solver_params_.end()) {
            in.get(params_[ParamKey::ORTHOGONALIZATION].name, orthogonalization_);
            solver_params->set(params_[ParamKey::ORTHOGONALIZATION].belos_name, orthogonalization_);
        }
        if (solver_params_.find(ParamKey::USE_SINGLE_RED) != solver_params_.end()) {
            in.get(params_[ParamKey::USE_SINGLE_RED].name, use_single_red_);
            solver_params->set(params_[ParamKey::USE_SINGLE_RED].belos_name, use_single_red_);
        }
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const std::shared_ptr<Preconditioner> &prec) {
        auto op = std::dynamic_pointer_cast<Operator<Vector>>(prec);

        if (!op) {
            assert(false && "IMPLEMENT ME");
        } else {
            auto tpetra_op = utopia::tpetra_operator(op);
            impl_->linear_problem->setLeftPrec(tpetra_op);
        }
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_programmatic_config() {
        using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

        auto solver_params = Teuchos::sublist(impl_->param_list, solver_type_, false);
        solver_params->set("Convergence Tolerance", static_cast<MagnitudeType>(this->rtol()));
        solver_params->set("Maximum Iterations", static_cast<int>(this->max_it()));
        solver_params->set("Output Style", 1);
        if (this->verbose()) {
            solver_params->set("Output Frequency", 1);
            solver_params->set("Verbosity",
                               Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary +
                                   Belos::StatusTestDetails);
        } else {
            solver_params->set("Verbosity", Belos::Errors + Belos::Warnings);
        }
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const std::string &pc_type,
                                                                   const std::string &pc_side) {
        if (pc_type.empty()) {
            enable_pc_ = false;
            enable_direct_pc_ = false;
        } else {
            std::string pc_type_lc(pc_type);
            std::transform(
                pc_type.begin(), pc_type.end(), pc_type_lc.begin(), [](unsigned char c) { return std::tolower(c); });
            if (pc_type_lc == "muelu") {
                enable_direct_pc_ = false;
            } else {
                direct_pc_type_ = pc_type;
                enable_direct_pc_ = true;
            }
            pc_side_ = pc_side;
            enable_pc_ = true;
        }
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_preconditioner(const std::shared_ptr<const Matrix> &op) {
        auto pc = Teuchos::RCP<typename Impl::OperatorType>();
        if (enable_direct_pc_) {
#ifdef UTOPIA_ENABLE_TRILINOS_IFPACK2
            const int solver_verbosity = impl_->param_list->get("Verbosity", Belos::Errors + Belos::Warnings);
            auto pc_params = Teuchos::sublist(impl_->param_list, direct_pc_type_, false);
            pc_params->set("Verbosity", solver_verbosity);
            Teuchos::RCP<typename Impl::IfPack2PrecType> direct_pc =
                Ifpack2::Factory::create<typename Impl::CrsMatrixType>(direct_pc_type_, raw_type(*op));
            direct_pc->setParameters(*pc_params);
            direct_pc->initialize();
            direct_pc->compute();
            pc = direct_pc;
#else   // UTOPIA_ENABLE_TRILINOS_IFPACK2
            m_utopia_error(
                "Cannot use Direct Preconditioner with BelosSolver, since Trilinos was built without Ifpack2 support");
#endif  // UTOPIA_ENABLE_TRILINOS_IFPACK2
        } else if (enable_pc_) {
#ifdef UTOPIA_ENABLE_TRILINOS_MUELU
            // Multigrid Hierarchy
            auto pc_params = Teuchos::sublist(impl_->param_list, "MueLu", false);
            pc_params->set("verbosity", this->verbose() ? "medium" : "none");
            Teuchos::RCP<typename Impl::MueLuPrecType> muelu_pc =
                MueLu::CreateTpetraPreconditioner((Teuchos::RCP<typename Impl::OperatorType>)raw_type(*op), *pc_params);

            pc = muelu_pc;
#else
            m_utopia_error(
                "Cannot use MueLu as preconditioner with BelosSolver, since Trilinos was built without MueLu support");
#endif  // UTOPIA_ENABLE_TRILINOS_MUELU
        }
        if (!pc.is_null()) {
            std::string pc_side_lc(pc_side_);
            std::transform(
                pc_side_.begin(), pc_side_.end(), pc_side_lc.begin(), [](unsigned char c) { return std::tolower(c); });
            if (pc_side_lc == "left") {
                impl_->linear_problem->setLeftPrec(pc);
            } else {
                impl_->linear_problem->setRightPrec(pc);
            }
        }
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::set_problem() {
        impl_->linear_problem->setProblem();

        auto solver_params = Teuchos::sublist(impl_->param_list, solver_type_, false);
        try {
            if (impl_->belos_solver.is_null()) {
                impl_->belos_solver = impl_->belos_factory.create(solver_type_, solver_params);
            } else {
                impl_->belos_solver->setParameters(solver_params);
            }
            impl_->belos_solver->setProblem(impl_->linear_problem);
        } catch (const std::exception &ex) {
            std::cerr << ex.what() << std::endl;
            assert(0);
        }

        if (this->verbose()) {
            std::cout << std::endl;
            std::cout << "Parameter list for '" << solver_type_ << "' solver:" << std::endl;
            impl_->belos_solver->getCurrentParameters()->print();
        }
    }

    template <typename Matrix, typename Vector>
    bool BelosSolver<Matrix, Vector, TRILINOS>::solve(const Operator<Vector> &op, const Vector &rhs, Vector &sol) {
        UTOPIA_TRACE_SCOPE("BelosSolver::solve(op)");

        if (empty(sol)) {
            sol.zeros(layout(rhs));
        }

        assert(!rhs.has_nan_or_inf());

        auto tpetra_op = tpetra_operator(utopia::make_ref(const_cast<Operator<Vector> &>(op)));

        impl_->linear_problem = Teuchos::rcp(new typename Impl::ProblemType(tpetra_op, raw_type(sol), raw_type(rhs)));

        set_programmatic_config();
        set_problem();
        // set_preconditioner(this->get_operator());

        assert(!(impl_->belos_solver.is_null()));
        impl_->belos_solver->solve();
        return true;
    }

    template <typename Matrix, typename Vector>
    void BelosSolver<Matrix, Vector, TRILINOS>::update(const Operator<Vector> &op) {
        // const auto layout_rhs = row_layout(op);

        // if (!initialized_ || !layout_rhs.same(layout_)) {
        //     init_memory(layout_rhs);
        // }
    }

}  // namespace utopia

#endif  // UTOPIA_BELOS_IMPL_HPP
#endif  // UTOPIA_ENABLE_TRILINOS_BELOS
