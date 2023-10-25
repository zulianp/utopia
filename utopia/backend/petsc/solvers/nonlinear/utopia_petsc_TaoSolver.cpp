#include "utopia_petsc_TaoSolver.hpp"
#include "utopia_Describable.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Types.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

#include <mpi.h>
#include <petsc/private/taoimpl.h>
#include "petsctao.h"

// Deprecated API is handled here
#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 17, 0)
#define Quircks_TaoSetHessian TaoSetHessianRoutine
#define Quircks_TaoSetGradient(input, b, c, d) TaoSetGradientRoutine(input, c, d)
#define Quircks_TaoSetObjective TaoSetObjectiveRoutine
#define Quircks_TaoSetSolution TaoSetInitialVector
#else
#define Quircks_TaoSetHessian TaoSetHessian
#define Quircks_TaoSetGradient(input, b, c, d) TaoSetGradient(input, b, c, d)
#define Quircks_TaoSetObjective TaoSetObjective
#define Quircks_TaoSetSolution TaoSetSolution
#endif

#define U_CHECKERR(ierr)               \
    {                                  \
        if ((ierr) != 0) return false; \
    }

namespace utopia {

    class TaoTypes : public Describable {
    public:
        inline static bool is_valid(const std::string &type, const bool verbose) {
            const auto &i = instance();
            bool valid = i.types_.find(type) != i.types_.end();

            if (!valid && verbose) {
                std::cerr << "Invalid tao type " << type << ". Valid types are: " << std::endl;
                i.describe(std::cerr);
            }

            return valid;
        }

        void describe(std::ostream &os) const override {
            for (const auto &t : types_) {
                os << t << " ";
            }

            os << std::endl;
        }

    private:
        std::set<std::string> types_;

        static inline const TaoTypes &instance() {
            static TaoTypes instance_;
            return instance_;
        }

        TaoTypes() {
            types_.insert(TAOLMVM);
            types_.insert(TAONLS);
            types_.insert(TAONTR);
            types_.insert(TAONTL);
            types_.insert(TAOCG);
            types_.insert(TAOTRON);
            types_.insert(TAOOWLQN);
            types_.insert(TAOBMRM);
            types_.insert(TAOBLMVM);

            types_.insert(TAOBQPIP);
            types_.insert(TAOGPCG);
            types_.insert(TAONM);
            types_.insert(TAOPOUNDERS);
            types_.insert(TAOLCL);
            types_.insert(TAOSSILS);
            types_.insert(TAOSSFLS);
            types_.insert(TAOASILS);
            types_.insert(TAOASFLS);
            types_.insert(TAOIPM);

            // TODO(zulianp): check if this is the right version
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 10, 3)

            types_.insert(TAOBQNLS);
            types_.insert(TAOBNCG);

            types_.insert(TAOBNLS);
            types_.insert(TAOBNTR);
            types_.insert(TAOBNTL);

            types_.insert(TAOBQNKLS);
            types_.insert(TAOBQNKTR);
            types_.insert(TAOBQNKTL);
#endif
        }
    };

    template <class Matrix, class Vector>
    static PetscErrorCode UtopiaTaoEvaluateObjective(Tao tao, Vec x, PetscReal *ret, void *ctx) {
        UTOPIA_UNUSED(tao);
        using Scalar = typename Function<Matrix, Vector>::Scalar;

        auto *fun = static_cast<Function<Matrix, Vector> *>(ctx);
        assert(fun);

        Vector utopia_x;
        convert(x, utopia_x);

        Scalar utopia_value = 0.;
        if (!fun->value(utopia_x, utopia_value)) {
            return 1;
        }

        *ret = utopia_value;
        return 0;
    }

    template <class Matrix, class Vector>
    static PetscErrorCode UtopiaTaoEvaluateGradient(Tao tao, Vec x, Vec g, void *ctx) {
        UTOPIA_UNUSED(tao);
        auto *fun = static_cast<Function<Matrix, Vector> *>(ctx);
        assert(fun);

        Vector utopia_x;
        convert(x, utopia_x);

        Vector utopia_g(layout(utopia_x));
        if (!fun->gradient(utopia_x, utopia_g)) {
            return 1;
        }

        convert(utopia_g, g);
        return 0;
    }

    template <class Matrix, class Vector>
    static PetscErrorCode UtopiaTaoFormHessian(Tao tao, Vec x, Mat H, Mat Hpre, void *ctx) {
        UTOPIA_UNUSED(tao);
        // PetscErrorCode ierr = 0;
        auto *fun = static_cast<Function<Matrix, Vector> *>(ctx);
        assert(fun);

        Vector utopia_x;

        Matrix utopia_H;
        Matrix utopia_Hpre;

        convert(x, utopia_x);
        utopia_H.wrap(H);
        utopia_Hpre.wrap(Hpre);

        if (!fun->hessian(utopia_x, utopia_H, utopia_Hpre)) {
            if (!fun->hessian(utopia_x, utopia_H)) {
                utopia_error("[Error] Failed to assemble Hessian.");
                assert(false);
                return 1;
            }

        } else {
            if (raw_type(utopia_Hpre) != Hpre) {
                // FIXME maybe add optimization options
                MatCopy(raw_type(utopia_Hpre), Hpre, DIFFERENT_NONZERO_PATTERN);
            }
        }

        if (raw_type(utopia_H) != H) {
            // FIXME maybe add optimization options
            MatCopy(raw_type(utopia_H), H, DIFFERENT_NONZERO_PATTERN);
        }

        return 0;
    }

    template <class Matrix, class Vector>
    bool UtopiaTaoSetUp(Tao tao, Function<Matrix, Vector> &fun) {
        fun.data()->init();
        if (!fun.initialize_hessian(*fun.data()->H, *fun.data()->H_pre)) {
            utopia_error("TaoSolver requires Function::initialize_hessian to be implemented.");
            assert(false);
            return false;
        }

        auto void_f = (static_cast<void *>(&fun));
        PetscErrorCode ierr = 0;
        if (fun.has_preconditioner()) {
            ierr = Quircks_TaoSetHessian(tao,
                                         raw_type(*fun.data()->H),
                                         raw_type(*fun.data()->H_pre),
                                         UtopiaTaoFormHessian<Matrix, Vector>,
                                         void_f);
            U_CHECKERR(ierr);

        } else {
            ierr = Quircks_TaoSetHessian(tao,
                                         raw_type(*fun.data()->H),
                                         raw_type(*fun.data()->H),
                                         (UtopiaTaoFormHessian<Matrix, Vector>),
                                         void_f);
            U_CHECKERR(ierr);
        }

        ierr = Quircks_TaoSetObjective(tao, UtopiaTaoEvaluateObjective<Matrix, Vector>, void_f);
        U_CHECKERR(ierr);

        ierr = Quircks_TaoSetGradient(tao, nullptr, (UtopiaTaoEvaluateGradient<Matrix, Vector>), void_f);
        U_CHECKERR(ierr);
        return false;
    }

    template <class Matrix, class Vector>
    class TaoSolver<Matrix, Vector>::Impl : public Configurable {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);
        using SizeType = UTOPIA_SIZE_TYPE(Vector);

        explicit Impl(MPI_Comm comm) : tao(nullptr), type_(TAOTRON) {
            assert(TaoTypes::is_valid(type_, true));
            init(comm);
        }

        Impl() : type_(TAOTRON) { assert(TaoTypes::is_valid(type_, true)); }

        ~Impl() override { destroy(); }

        void init(MPI_Comm comm) {
            destroy();
            assert(TaoTypes::is_valid(type_, true));

            auto ierr = TaoCreate(comm, &tao);
            assert(ierr == 0);
            ierr = TaoSetType(tao, type_.c_str());
            assert(ierr == 0);

            UTOPIA_UNUSED(ierr);
        }

        void set_from_options() {
            assert(initialized());
            auto ierr = TaoSetFromOptions(tao);
            assert(ierr == 0);
            UTOPIA_UNUSED(ierr);
        }

        void destroy() {
            if (tao) {
                TaoDestroy(&tao);
                tao = nullptr;
            }
        }

        inline std::string get_type() const {
            if (initialized()) {
                // TODO(zulianp): check if this is the right version
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 10, 2)
                TaoType type;
#else
                const TaoType type;
#endif
                TaoGetType(tao, &type);

                assert(TaoTypes::is_valid(type, true));
                return type;
            }

            assert(TaoTypes::is_valid(type_, true));

            return type_;
        }

        void set_type(const std::string &type) {
            assert(TaoTypes::is_valid(type, true));

            type_ = type;
            if (tao) {
                TaoSetType(tao, type_.c_str());
            }
        }

        void set_bounds(const Vector &lb, const Vector &ub) {
            assert(initialized());

            if (tao == nullptr) {
                std::cerr << "[Error] attempt to set bounds to uninitialized tao" << std::endl;
                return;
            }

            auto ierr = TaoSetVariableBounds(tao, raw_type(lb), raw_type(ub));
            assert(ierr == 0);
            UTOPIA_UNUSED(ierr);
        }

        bool get_ksp(KSP *ksp) {
            assert(initialized());

            if (tao == nullptr) {
                std::cerr << "[Error] attempt to set ksp to uninitialized tao" << std::endl;
                return false;
            }

            PetscErrorCode ierr = TaoGetKSP(tao, ksp);
            assert(ierr == 0);
            return ierr == 0;
        }

        void set_tol(const Scalar gatol, const Scalar grtol, const Scalar gttol, const SizeType maxits) {
            assert(initialized());

            if (tao == nullptr) {
                std::cerr << "[Error] attempt to set tol to uninitialized tao" << std::endl;
                return;
            }

            auto ierr = TaoSetTolerances(tao, gatol, grtol, gttol);
            assert(ierr == 0);
            ierr = TaoSetMaximumIterations(tao, maxits);
            assert(ierr == 0);
            UTOPIA_UNUSED(ierr);
        }

        void set_monitor() {
            const char monfilename[7] = "stdout";
            PetscViewer monviewer;
            PetscViewerASCIIOpen(communicator(), monfilename, &monviewer);
            TaoSetMonitor(
                tao, TaoDefaultSMonitor, monviewer, reinterpret_cast<PetscErrorCode (*)(void **)>(PetscViewerDestroy));
        }

        void read(Input &in) override {
            std::string type;
            in.get("type", type);
            in.get("tao_type", type);

            if (!type.empty() && TaoTypes::is_valid(type, true)) {
                set_type(type);
            }
        }

        inline bool initialized(const SizeType &n_global) const {
            if (tao != nullptr) {
                PetscInt size;
                VecGetSize((*tao).solution, &size);
                return (size == n_global);
            }
            { return false; }
        }

        inline bool initialized() const { return (tao != nullptr); }

        void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &solver) {
            assert(tao);

            KSP ksp;
            get_ksp(&ksp);

            if (ksp != nullptr) {
                build_ksp(solver, ksp);
            }
        }

        void set_function(Function<Matrix, Vector> &fun) { UtopiaTaoSetUp(tao, fun); }

        inline bool solve(Vector &x) {
            PetscErrorCode ierr = 0;
            Quircks_TaoSetSolution(tao, raw_type(x));
            ierr = TaoSolve(tao);
            U_CHECKERR(ierr);

            PetscInt iterate;
            PetscReal f;
            PetscReal gnorm;
            PetscReal cnorm;
            PetscReal xdiff;
            TaoConvergedReason reason;
            TaoGetSolutionStatus(tao, &iterate, &f, &gnorm, &cnorm, &xdiff, &reason);

            // if(this->verbose()) {
            // utopia::out() <<"iterate: " << iterate << std::endl;
            // utopia::out() <<"f: " << f << std::endl;
            // utopia::out() <<"gnorm: " << gnorm << std::endl;
            // utopia::out() <<"cnorm: " << cnorm << std::endl;
            // utopia::out() <<"xdiff: " << xdiff << std::endl;
            // utopia::out() <<"reason: " << reason << std::endl;
            // }

            if (reason < 0) {
                utopia_warning("> Failed to converge");

                std::cerr << "gnorm: " << gnorm << std::endl;
                std::cerr << "cnorm: " << cnorm << std::endl;
                std::cerr << "xdiff: " << xdiff << std::endl;
            }

            return reason >= 0;
        }

        inline void get_sol_status(PetscInt &iterates, TaoConvergedReason &reason) {
            PetscReal f;
            PetscReal gnorm;
            PetscReal cnorm;
            PetscReal xdiff;
            TaoGetSolutionStatus(tao, &iterates, &f, &gnorm, &cnorm, &xdiff, &reason);
        }

        inline bool smooth(Vector &x) {
            PetscErrorCode ierr = 0;
            Quircks_TaoSetSolution(tao, raw_type(x));
            ierr = TaoSolve(tao);
            U_CHECKERR(ierr);
            return true;
        }

        inline MPI_Comm communicator() const {
            assert(initialized());
            MPI_Comm comm = PetscObjectComm((PetscObject)tao);
            assert(comm != MPI_COMM_NULL);
            return comm;
        }

    private:
        Tao tao{nullptr};
        std::string type_;
    };

    template <class Matrix, class Vector>
    TaoSolver<Matrix, Vector>::TaoSolver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &linear_solver)
        : NewtonBase<Matrix, Vector>(linear_solver) {
        impl_ = utopia::make_unique<Impl>();
    }

    template <class Matrix, class Vector>
    TaoSolver<Matrix, Vector>::TaoSolver() : NewtonBase<Matrix, Vector>(nullptr) {
        impl_ = utopia::make_unique<Impl>();
    }

    template <class Matrix, class Vector>
    TaoSolver<Matrix, Vector>::~TaoSolver() = default;

    template <class Matrix, class Vector>
    void TaoSolver<Matrix, Vector>::set_type(const std::string &type) {
        impl_->set_type(type);
    }

    template <class Matrix, class Vector>
    void TaoSolver<Matrix, Vector>::read(Input &in) {
        in.get("linear_solver", [&](Input &) {
            auto ls = std::make_shared<OmniLinearSolver<Matrix, Vector>>();
            this->set_linear_solver(ls);
        });

        NewtonBase<Matrix, Vector>::read(in);
        // VariableBoundSolverInterface<Vector>::read(in);
        impl_->read(in);
    }

    template <class Matrix, class Vector>
    void TaoSolver<Matrix, Vector>::print_usage(std::ostream &os) const {
        NewtonBase<Matrix, Vector>::print_usage(os);
        this->print_param_usage(os, "type", "string", "Type of tao solver.", "-");
    }

    template <class Matrix, class Vector>
    bool TaoSolver<Matrix, Vector>::solve(Function<Matrix, Vector> &fun, Vector &x) {
        UTOPIA_TRACE_SCOPE("TaoSolver::solve");

        init(fun, x);
        this->init_solver("Tao Solver", {""});
        auto flg = impl_->solve(x);

        PetscInt iterates;
        TaoConvergedReason reason;

        impl_->get_sol_status(iterates, reason);
        this->exit_solver(iterates, reason);
        this->print_statistics(iterates);

        return flg;
    }

    template <class Matrix, class Vector>
    bool TaoSolver<Matrix, Vector>::smooth(Function<Matrix, Vector> &fun, Vector &x) {
        init(fun, x);
        return impl_->smooth(x);
    }

    template <class Matrix, class Vector>
    void TaoSolver<Matrix, Vector>::init(Function<Matrix, Vector> &fun, Vector &x) {
        UTOPIA_TRACE_SCOPE("TaoSolver::init");

        if (!impl_->initialized(size(x))) {
            impl_->init(x.comm().get());

            if (this->linear_solver()) {
                impl_->set_linear_solver(this->linear_solver());
            }
        }

        impl_->set_tol(this->atol(), this->rtol(), this->stol(), this->max_it());

        if (this->has_bound()) {
            this->fill_empty_bounds(layout(x));
            impl_->set_bounds(this->get_lower_bound(), this->get_upper_bound());
        }

        impl_->set_function(fun);

        if (this->verbose()) {
            impl_->set_monitor();
        }

        impl_->set_from_options();
    }

    template <class Matrix, class Vector>
    bool TaoSolver<Matrix, Vector>::get_ksp(KSP *ksp) {
        return impl_->get_ksp(ksp);
    }

    template <class Matrix, class Vector>
    TaoSolver<Matrix, Vector> *TaoSolver<Matrix, Vector>::clone() const {
        // FIXME make proper deep clone
        auto cloned = utopia::make_unique<TaoSolver>();

        cloned->set_type(cloned->impl_->get_type());

        if (this->impl_->initialized()) {
            cloned->impl_->init(impl_->communicator());
        }

        return cloned.release();
    }

    template class TaoSolver<PetscMatrix, PetscVector>;
    // FIXME
    // template class TaoSolver<PetscMatrix, PetscVector>;
}  // namespace utopia

#undef U_CHECKERR
