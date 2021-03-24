#ifndef UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP
#define UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP

#include "utopia_petsc_KSPSolver.hpp"

#include "utopia_Core.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_SolutionStatus.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_petsc_KSPTypes.hpp"

#include <petsc/private/kspimpl.h>
#include <petscksp.h>
#include <petscpc.h>
#include <petscsys.h>
#include <algorithm>

namespace utopia {

    class KSPLog {
    public:
        Vec x_k_1{nullptr};
        Vec x_k_2{nullptr};

        KSPLog() = default;

        inline bool initialized() const { return x_k_1 != nullptr; }

        void init_from(Vec x) {
            destroy();

            VecDuplicate(x, &x_k_1);
            VecDuplicate(x, &x_k_2);
        }

        ~KSPLog() { destroy(); }

        inline void destroy() {
            if (x_k_1) {
                VecDestroy(&x_k_1);
                x_k_1 = nullptr;
            }

            if (x_k_2) {
                VecDestroy(&x_k_2);
                x_k_2 = nullptr;
            }
        }
    };

    template <typename Matrix, typename Vector>
    class KSPSolver<Matrix, Vector, PETSC>::Impl {
    public:
        Impl(MPI_Comm comm) : ksp_(nullptr), owner_(true) {
            init(comm);

            ksp_type(KSPTypes::instance().ksp(0));
            pc_type(KSPTypes::instance().pc(0));
            solver_package(KSPTypes::instance().package(0));
        }

        Impl(KSP &ksp, const bool owner = false) : ksp_(ksp), owner_(owner) { set_command_line_parameters(); }

        /* @brief      Sets the choice of direct solver.
         *             Please note, in petsc, direct solver is used as preconditioner alone, with proper settings.
         *
         * @param[in]  PCType  The type of direct solver.
         */
        inline void pc_type(const std::string &pc_type) {
            if (KSPTypes::instance().is_pc_valid(pc_type)) {
                PC pc;
                KSPGetPC(ksp_, &pc);
                PCSetType(pc, pc_type.c_str());
            }
        }

        /**
         * @brief      Sets KSP type
         */
        inline void ksp_type(const std::string &ksp_type) {
            if (KSPTypes::instance().is_ksp_valid(ksp_type)) {
                KSPSetType(ksp_, ksp_type.c_str());
            }
        }

        /**
         * @brief      Sets solver package for choice of direct solver.
         *
         * @param[in]  SolverPackage  The solver package.
         */
        inline void solver_package(const std::string &package) {
            if (KSPTypes::instance().is_solver_package_valid(package)) {
                PetscErrorCode ierr;
                UTOPIA_UNUSED(ierr);
                PC pc;
                ierr = KSPGetPC(ksp_, &pc);
                assert(ierr == 0);

#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 9, 0)
                ierr = PCFactorSetMatSolverPackage(pc, package.c_str());
                assert(ierr == 0);
#else
                ierr = PCFactorSetMatSolverType(pc, package.c_str());
                assert(ierr == 0);
#endif
            }
            // else {
            //     std::cerr << "Invalid solver package \"" << package << "\"" << std::endl;
            // }
        }

        /**
         * @brief      Returns type of direct solver.
         */
        inline PCType pc_type() const {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            PC pc;
            PCType ret;
            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);
            ierr = PCGetType(pc, &ret);
            assert(ierr == 0);
            return ret;
        }

        inline bool has_shell_pc() const { return std::string(PCSHELL) == pc_type(); }

        /**
         * @brief      Returns ksp package type
         */
        inline KSPType ksp_type() const {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            KSPType ret;
            ierr = KSPGetType(ksp_, &ret);
            assert(ierr == 0);
            return ret;
        }

        /**
         * @brief      Returns type of solver package.
         */
        inline std::string solver_package() const {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);

            PC pc;
            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);

#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 9, 0)
            const MatSolverPackage stype;
            ierr = PCFactorGetMatSolverPackage(pc, &stype);
            assert(ierr == 0);
            if (stype)
                return stype;
            else
                return "";
#else
            MatSolverType stype;
            ierr = PCFactorGetMatSolverType(pc, &stype);
            assert(ierr == 0);
            if (stype)
                return stype;
            else
                return "";
#endif

            // static const char * ret = " ";
            // return ret;
        }

        inline bool is_null() const { return ksp_ == nullptr; }

        ~Impl() { destroy(); }

        inline void destroy() {
            if (ksp_) {
                if (owner_) {
                    KSPDestroy(&ksp_);
                }

                ksp_ = nullptr;
            }
        }

        inline MPI_Comm communicator() const {
            MPI_Comm comm = PetscObjectComm((PetscObject)ksp_);
            assert(comm != MPI_COMM_NULL);
            return comm;
        }

        void describe() const {
            PetscBool bool_value;
            KSPGetInitialGuessNonzero(ksp_, &bool_value);
            utopia::out() << "-------------------------------------------" << std::endl;
            utopia::out() << "KSPGetInitialGuessNonzero: " << bool_value << std::endl;
            utopia::out() << "PCType:                    " << pc_type() << std::endl;
            utopia::out() << "KSPType:                   " << ksp_type() << std::endl;
            utopia::out() << "-------------------------------------------" << std::endl;

            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            ierr = KSPView(ksp_, PETSC_VIEWER_STDOUT_(communicator()));
            assert(ierr == 0);
        }

        inline void set_from_options() { KSPSetFromOptions(ksp_); }

        inline void init(const MPI_Comm comm) {
            destroy();
            KSPCreate(comm, &ksp_);
            KSPSetFromOptions(ksp_);
            KSPSetComputeSingularValues(ksp_, PETSC_FALSE);
            // KSPSetNormType(ksp_, KSP_NORM_UNPRECONDITIONED);
        }

        inline void set_initial_guess_non_zero(const bool val) {
            if (val) {
                KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
            } else {
                KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);
            }
        }

        void pc_copy_settings(PC &other_pc) const {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            PC this_pc;
            PCType type;

            ierr = KSPGetPC(ksp_, &this_pc);
            assert(ierr == 0);

            ierr = PCGetType(this_pc, &type);
            assert(ierr == 0);
            ierr = PCSetType(other_pc, type);
            assert(ierr == 0);

#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 9, 0)
            const MatSolverPackage stype;
            ierr = PCFactorGetMatSolverPackage(this_pc, &stype);
            assert(ierr == 0);
            ierr = PCFactorSetMatSolverPackage(other_pc, stype);
            assert(ierr == 0);
#else
            MatSolverType stype;
            ierr = PCFactorGetMatSolverType(this_pc, &stype);
            assert(ierr == 0);
            ierr = PCFactorSetMatSolverType(other_pc, stype);
            assert(ierr == 0);
#endif

            PetscBool flg_is_redundant;
            PetscObjectTypeCompare((PetscObject)other_pc, PCREDUNDANT, &flg_is_redundant);

            if (flg_is_redundant) {
                // there is no function to get number, so it can not be coppied....
                // PCRedundantSetNumber(other_pc, number);

                // let us copy at least ksp and pc types

                // setting up inner solver
                KSP inner_ksp_other, inner_ksp_this;
                PCRedundantGetKSP(other_pc, &inner_ksp_other);
                PCRedundantGetKSP(this_pc, &inner_ksp_this);

                KSPType inner_ksp_type;
                KSPGetType(inner_ksp_this, &inner_ksp_type);
                KSPSetType(inner_ksp_other, inner_ksp_type);

                PC innner_PC_this, innner_PC_other;
                KSPGetPC(inner_ksp_this, &innner_PC_this);
                KSPGetPC(inner_ksp_other, &innner_PC_other);

                PCType inner_pc_type;
                PCGetType(innner_PC_this, &inner_pc_type);
                PCSetType(innner_PC_other, inner_pc_type);
            }
        }

        void copy_settings_from(const Impl &other) { other.copy_settings_to(ksp_); }

        void copy_settings_to(KSP &other_ksp) const {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            PetscBool bool_value;
            KSPType type;
            PetscReal rtol;
            PetscReal abstol;
            PetscReal dtol;
            PetscInt maxits;
            PC other_pc;

            KSPSetFromOptions(other_ksp);

            ierr = KSPGetComputeSingularValues(ksp_, &bool_value);
            assert(ierr == 0);
            ierr = KSPSetComputeSingularValues(other_ksp, bool_value);
            assert(ierr == 0);

            ierr = KSPGetType(ksp_, &type);
            assert(ierr == 0);
            ierr = KSPSetType(other_ksp, type);
            assert(ierr == 0);

            ierr = KSPGetInitialGuessNonzero(ksp_, &bool_value);
            assert(ierr == 0);
            ierr = KSPSetInitialGuessNonzero(other_ksp, bool_value);
            assert(ierr == 0);

            ierr = KSPGetTolerances(ksp_, &rtol, &abstol, &dtol, &maxits);
            assert(ierr == 0);
            ierr = KSPSetTolerances(other_ksp, rtol, abstol, dtol, maxits);
            assert(ierr == 0);

            ierr = KSPGetPC(other_ksp, &other_pc);

            pc_copy_settings(other_pc);
        }

        void attach_shell_preconditioner(PetscErrorCode (*apply)(PC, Vec, Vec),
                                         void *ctx,
                                         PetscErrorCode (*setup)(PC),
                                         PetscErrorCode (*destroy)(PC)) {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            PC pc;

            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);
            ierr = PCSetType(pc, PCSHELL);
            assert(ierr == 0);
            ierr = PCShellSetApply(pc, apply);
            assert(ierr == 0);

            ierr = PCShellSetName(pc, "Utopia Preconditioner");
            assert(ierr == 0);

            if (setup) {
                ierr = PCShellSetSetUp(pc, setup);
                assert(ierr == 0);
            }

            if (destroy) {
                ierr = PCShellSetSetUp(pc, destroy);
                assert(ierr == 0);
            }

            if (ctx) {
                ierr = PCShellSetContext(pc, ctx);
                assert(ierr == 0);
            }

            // usefull methods
            // PCGetOperatorsSet(PC pc,PetscBool  *mat,PetscBool  *pmat)
            // PetscErrorCode  PCGetOperators(PC pc,Mat *Amat,Mat *Pmat)
            // PetscErrorCode  PCShellSetApplyRichardson(PC pc,PetscErrorCode
            // (*apply)(PC,Vec,Vec,Vec,PetscReal,PetscReal,PetscReal,PetscInt,PetscBool,PetscInt*,PCRichardsonConvergedReason*))
            // PCGetReusePreconditioner(PC pc,PetscBool *flag)
        }

        void update(const Matrix &mat) {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);

            ierr = KSPSetOperators(ksp_, raw_type(mat), raw_type(mat));
            assert(ierr == 0);

            // should not be necessary
            // ierr = KSPSetUp(ksp_);                                      assert(ierr == 0);
        }

        void update(const Matrix &mat, const Matrix &prec) {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);

            ierr = KSPSetOperators(ksp_, raw_type(mat), raw_type(prec));
            assert(ierr == 0);

            // should not be necessary
            // ierr = KSPSetUp(ksp_);                                       assert(ierr == 0);
        }

        bool smooth(const SizeType sweeps, const Vector &rhs, Vector &x) {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            KSPNormType normtype;
            PetscReal rtol;
            PetscReal abstol;
            PetscReal dtol;
            PetscInt maxits;
            void *ctx;

            // back-up settings
            ierr = KSPGetNormType(ksp_, &normtype);
            assert(ierr == 0);
            ierr = KSPGetTolerances(ksp_, &rtol, &abstol, &dtol, &maxits);
            assert(ierr == 0);

            // set up smoothing things
            ierr = KSPSetTolerances(ksp_, 0., 0., PETSC_DEFAULT, sweeps);
            assert(ierr == 0);
            ierr = KSPSetNormType(ksp_, KSP_NORM_NONE);
            assert(ierr == 0);
            ierr = KSPSetConvergenceTest(ksp_, KSPConvergedSkip, nullptr, nullptr);
            assert(ierr == 0);

            // perform smoothing
            ierr = KSPSolve(ksp_, raw_type(rhs), raw_type(x));
            assert(ierr == 0);

            // reset-solver-stuff
            ierr = KSPSetNormType(ksp_, normtype);
            assert(ierr == 0);
            ierr = KSPSetTolerances(ksp_, rtol, abstol, dtol, maxits);
            assert(ierr == 0);
            ierr = KSPConvergedDefaultCreate(&ctx);
            assert(ierr == 0);
            ierr = KSPSetConvergenceTest(ksp_, KSPConvergedDefault, ctx, KSPConvergedDefaultDestroy);
            assert(ierr == 0);
            return true;
        }

        void set_tolerances(const PetscReal rtol, const PetscReal atol, const PetscReal dtol, const PetscInt max_it) {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            ierr = KSPSetTolerances(ksp_, rtol, atol, dtol, max_it);
            assert(ierr == 0);
        }

        void number_of_subdomains(const PetscInt &n_blocks) {
            PetscErrorCode ierr;
            PC pc;
            PCType pc_type;

            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);
            ierr = PCGetType(pc, &pc_type);
            assert(ierr == 0);

            PetscBool flg_bjacobi = PETSC_FALSE, flg_asm = PETSC_FALSE;
            PetscObjectTypeCompare((PetscObject)pc, PCBJACOBI, &flg_bjacobi);

            if (!flg_bjacobi) {
                PetscObjectTypeCompare((PetscObject)pc, PCASM, &flg_asm);
            }

            if (flg_bjacobi) {
                PCBJacobiSetTotalBlocks(pc, n_blocks, nullptr);
            } else if (flg_asm) {
                PCGASMSetTotalSubdomains(pc, n_blocks);
            } else {
                utopia_error("Number of subdomain can be set only for PCBJACOBI or PCASM.");
            }

            UTOPIA_UNUSED(ierr);
        }

        void norm_type(const std::string &norm_type) {
            if (norm_type == "preconditioned") {
                KSPSetNormType(ksp_, KSP_NORM_PRECONDITIONED);
            } else if (norm_type == "unpreconditioned") {
                KSPSetNormType(ksp_, KSP_NORM_UNPRECONDITIONED);
            } else if (norm_type == "none") {
                KSPSetNormType(ksp_, KSP_NORM_NONE);
            } else if (norm_type == "natural") {
                KSPSetNormType(ksp_, KSP_NORM_NATURAL);
            } else {
                utopia_warning("KSP::norm_type:: norm not supported... \n");
            }
        }

        void overlap(const PetscInt &n_overlap) {
            PetscErrorCode ierr;
            PC pc;
            PCType pc_type;

            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);
            ierr = PCGetType(pc, &pc_type);
            assert(ierr == 0);

            PetscBool flg_asm = PETSC_FALSE;
            PetscObjectTypeCompare((PetscObject)pc, PCASM, &flg_asm);

            if (flg_asm) {
                PCGASMSetOverlap(pc, n_overlap);
            } else {
                utopia_error("Overlap can be set only for PCASM");
                return;
            }

            UTOPIA_UNUSED(ierr);
        }

        void sub_ksp_pc_type(const std::string &sub_ksp_type, const std::string &sub_pc_type) {
            if (!ksp_->setupstage) {
                // ksp setup has to be called before setting up preconditioner
                KSPSetUp(ksp_);
            }

            PetscErrorCode ierr;
            PC pc;
            PCType pc_type;

            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);
            ierr = PCGetType(pc, &pc_type);
            assert(ierr == 0);

            PetscBool flg_bjacobi = PETSC_FALSE, flg_asm = PETSC_FALSE;
            PetscObjectTypeCompare((PetscObject)pc, PCBJACOBI, &flg_bjacobi);

            PetscInt first, nlocal;
            KSP *subksp;

            if (!flg_bjacobi) {
                PetscObjectTypeCompare((PetscObject)pc, PCASM, &flg_asm);
            }

            if (flg_bjacobi) {
                PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
            } else if (flg_asm) {
                PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
            } else {
                utopia_error("Sub_ksp_pc_type can be set only for PCBJACOBI and PCASM");
                return;
            }

            for (auto i = 0; i < nlocal; i++) {
                KSPSetType(subksp[i], sub_ksp_type.c_str());
                PC subpc;
                KSPGetPC(subksp[i], &subpc);
                PCSetType(subpc, sub_pc_type.c_str());
            }

            UTOPIA_UNUSED(ierr);
        }

        void sub_ksp_type(const std::string &sub_ksp_type) {
            if (!ksp_->setupstage) {
                utopia_error("sub_ksp_type can be only called after update(). ");
                return;
            }

            PetscErrorCode ierr;
            PC pc;
            PCType pc_type;

            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);
            ierr = PCGetType(pc, &pc_type);
            assert(ierr == 0);

            PetscBool flg_bjacobi = PETSC_FALSE, flg_asm = PETSC_FALSE;
            PetscObjectTypeCompare((PetscObject)pc, PCBJACOBI, &flg_bjacobi);

            PetscInt first, nlocal;
            KSP *subksp;

            if (!flg_bjacobi) {
                PetscObjectTypeCompare((PetscObject)pc, PCASM, &flg_asm);
            }

            if (flg_bjacobi) {
                PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
            } else if (flg_asm) {
                PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
            } else {
                utopia_error("Sub_ksp_pc_type can be set only for PCBJACOBI and PCASM");
                return;
            }

            for (auto i = 0; i < nlocal; i++) {
                KSPSetType(subksp[i], sub_ksp_type.c_str());
            }

            UTOPIA_UNUSED(ierr);
        }

        void sub_pc_type(const std::string &sub_pc_type) {
            if (!ksp_->setupstage) {
                utopia_error("sub_ksp_pc_type can be only called after update(). ");
                return;
            }

            PetscErrorCode ierr;
            PC pc;
            PCType pc_type;

            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);
            ierr = PCGetType(pc, &pc_type);
            assert(ierr == 0);

            PetscBool flg_bjacobi = PETSC_FALSE, flg_asm = PETSC_FALSE;
            PetscObjectTypeCompare((PetscObject)pc, PCBJACOBI, &flg_bjacobi);

            PetscInt first, nlocal;
            KSP *subksp;

            if (!flg_bjacobi) {
                PetscObjectTypeCompare((PetscObject)pc, PCASM, &flg_asm);
            }

            if (flg_bjacobi) {
                PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
            } else if (flg_asm) {
                PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
            } else {
                utopia_error("Sub_ksp_pc_type can be set only for PCBJACOBI and PCASM");
                return;
            }

            for (auto i = 0; i < nlocal; i++) {
                PC subpc;
                KSPGetPC(subksp[i], &subpc);
                PCSetType(subpc, sub_pc_type.c_str());
            }

            UTOPIA_UNUSED(ierr);
        }

        void sub_solver_package(const std::string &sub_pc_type) {
            if (!ksp_->setupstage) {
                utopia_error("sub_ksp_pc_type can be only called after update(). ");
                return;
            }

            PetscErrorCode ierr;
            PC pc;
            PCType pc_type;

            ierr = KSPGetPC(ksp_, &pc);
            assert(ierr == 0);
            ierr = PCGetType(pc, &pc_type);
            assert(ierr == 0);

            PetscBool flg_bjacobi = PETSC_FALSE, flg_asm = PETSC_FALSE;
            PetscObjectTypeCompare((PetscObject)pc, PCBJACOBI, &flg_bjacobi);

            PetscInt first, nlocal;
            KSP *subksp;

            if (!flg_bjacobi) {
                PetscObjectTypeCompare((PetscObject)pc, PCASM, &flg_asm);
            }

            if (flg_bjacobi) {
                PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
            } else if (flg_asm) {
                PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
            } else {
                utopia_error("Sub_ksp_pc_type can be set only for PCBJACOBI and PCASM");
                return;
            }

            for (auto i = 0; i < nlocal; i++) {
                PC subpc;
                KSPGetPC(subksp[i], &subpc);
#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 9, 0)
                PCFactorSetMatSolverPackage(subpc, sub_pc_type.c_str());
#else
                PCFactorSetMatSolverType(subpc, sub_pc_type.c_str());
#endif
            }

            UTOPIA_UNUSED(ierr);
        }

        void set_monitor(PetscErrorCode (*monitor)(KSP, PetscInt, PetscReal, void *),
                         void *mctx,
                         PetscErrorCode (*monitordestroy)(void **)) {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            ierr = KSPMonitorSet(ksp_, monitor, mctx, monitordestroy);
            assert(ierr == 0);
            UTOPIA_UNUSED(ierr);
        }

        bool apply(const Vector &b, Vector &x) {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);

            if (empty(x) || (size(b).get(0) != size(x).get(0))) {
                x.zeros(layout(b));
            }

            assert(local_size(b).get(0) == local_size(x).get(0));

            KSPConvergedReason reason;
            ierr = KSPSolve(ksp_, raw_type(b), raw_type(x));
            assert(ierr == 0);
            ierr = KSPGetConvergedReason(ksp_, &reason);
            assert(ierr == 0);

            // if(reason < 0) {

            //     utopia_warning(
            //         "ksp apply returned " + std::to_string(reason) + " = " + converged_str(reason) +
            //         " ksp_type=" + ksp_type() + " pc_type=" + pc_type() + " solver_package=" + solver_package());
            // }

            return reason >= 0;
        }

        void solution_status(SolutionStatus &status) {
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            PetscInt its;
            KSPConvergedReason reason;
            PetscReal rnorm;

            ierr = KSPGetIterationNumber(ksp_, &its);
            assert(ierr == 0);
            ierr = KSPGetConvergedReason(ksp_, &reason);
            assert(ierr == 0);
            ierr = KSPGetResidualNorm(ksp_, &rnorm);
            assert(ierr == 0);

            status.iterates = its;
            status.reason = reason;
            status.gradient_norm = rnorm;
        }

        void set_command_line_parameters() {
            // petsc command line options
            char name_[1024];
            PetscBool flg;

#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 7, 0)
            PetscOptionsGetString(nullptr, "-ksp_type", name_, 1024, &flg);
            if (flg) {
                ksp_type(name_);
            }

            PetscOptionsGetString(nullptr, "-pc_type", name_, 1024, &flg);
            if (flg) {
                pc_type(name_);
            }

            PetscOptionsGetString(nullptr, "-pc_factor_mat_solver_package", name_, 1024, &flg);
            if (flg) {
                solver_package(name_);
            }

#else
            PetscOptionsGetString(nullptr, nullptr, "-ksp_type", name_, 1024, &flg);
            if (flg) {
                ksp_type(name_);
            }

            PetscOptionsGetString(nullptr, nullptr, "-pc_type", name_, 1024, &flg);
            if (flg) {
                pc_type(name_);
            }

            PetscOptionsGetString(nullptr, nullptr, "-pc_factor_mat_solver_package", name_, 1024, &flg);
            if (flg) {
                solver_package(name_);
            }

            PetscOptionsGetString(nullptr, nullptr, "-pc_factor_mat_solver_type", name_, 1024, &flg);
            if (flg) {
                solver_package(name_);
            }
#endif
        }

        KSP &implementation() { return ksp_; }

        const KSP &implementation() const { return ksp_; }

    private:
        KSP ksp_;
        bool owner_;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // FIXME defaults to PETSC_COMM_WORLD
    template <typename Matrix, typename Vector>
    KSPSolver<Matrix, Vector, PETSC>::KSPSolver() : ksp_(utopia::make_unique<Impl>(PETSC_COMM_WORLD)) {
        ksp_type("bcgs");
        pc_type("jacobi");
        ksp_->set_from_options();
        ksp_->set_initial_guess_non_zero(true);
    }

    template <typename Matrix, typename Vector>
    KSPSolver<Matrix, Vector, PETSC>::KSPSolver(
        std::unique_ptr<typename utopia::KSPSolver<Matrix, Vector, PETSC>::Impl> &&w)
        : ksp_(std::move(w)) {}

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::wrap(KSP &ksp) {
        ksp_ = utopia::make_unique<Impl>(ksp, false);
    }

    template <typename Matrix, typename Vector>
    KSPSolver<Matrix, Vector, PETSC>::~KSPSolver() = default;

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::set_initial_guess_non_zero(const bool val) {
        ksp().set_initial_guess_non_zero(val);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::factor_set_pivot_in_blocks(const bool val) {
        PC pc;
        auto ierr = KSPGetPC(ksp_->implementation(), &pc);
        assert(ierr == 0);
        UTOPIA_UNUSED(ierr);

        PCFactorSetPivotInBlocks(pc, val ? PETSC_TRUE : PETSC_FALSE);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::pc_type(const std::string &pc_type) {
        ksp_->pc_type(pc_type);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::ksp_type(const std::string &ksp_type) {
        ksp_->ksp_type(ksp_type);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::solver_package(const std::string &package) {
        ksp_->solver_package(package);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::number_of_subdomains(const SizeType &n) {
        return this->ksp_->number_of_subdomains(n);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::norm_type(const std::string &norm_type) {
        return this->ksp_->norm_type(norm_type);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::overlap(const SizeType &n) {
        return this->ksp_->overlap(n);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::sub_ksp_pc_type(const std::string ksp_type, const std::string pc_type) {
        return this->ksp_->sub_ksp_pc_type(ksp_type, pc_type);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::sub_ksp_type(const std::string type) {
        return this->ksp_->sub_ksp_type(type);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::sub_pc_type(const std::string type) {
        return this->ksp_->sub_pc_type(type);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::sub_solver_package(const std::string type) {
        return this->ksp_->sub_solver_package(type);
    }

    template <typename Matrix, typename Vector>
    std::string KSPSolver<Matrix, Vector, PETSC>::pc_type() const {
        return this->ksp_->pc_type();
    }

    template <typename Matrix, typename Vector>
    std::string KSPSolver<Matrix, Vector, PETSC>::ksp_type() const {
        return this->ksp_->ksp_type();
    }

    template <typename Matrix, typename Vector>
    std::string KSPSolver<Matrix, Vector, PETSC>::solver_package() const {
        return this->ksp_->solver_package();
    }

    template <typename Matrix, typename Vector>
    bool KSPSolver<Matrix, Vector, PETSC>::apply(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("KSPSolver::apply");

        ksp_->set_tolerances(this->rtol(), this->atol(), PETSC_DEFAULT, this->max_it());

        bool flg = ksp_->apply(b, x);

        ksp_->solution_status(this->solution_status_);

        // is this proper place to do so???
        // this->set_ksp_options(ksp_->implementation());

        UTOPIA_TRACE_REGION_END("KSPSolver::apply");
        return flg;
    }

    template <typename Matrix, typename Vector>
    bool KSPSolver<Matrix, Vector, PETSC>::smooth(const Vector &rhs, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("KSPSolver::smooth");

        bool ok = ksp_->smooth(this->sweeps(), rhs, x);

        UTOPIA_TRACE_REGION_END("KSPSolver::smooth");
        return ok;
    }

    template <typename Matrix, typename Vector>
    bool KSPSolver<Matrix, Vector, PETSC>::must_compute_cond_number() const {
        static const std::vector<std::string> types = {"bcgs", "cg", "gmres"};
        return (this->verbose() && in_array(ksp_->ksp_type(), types));
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::set_ksp_options(KSP &ksp) {
        ksp_->copy_settings_to(ksp);
        set_monitor_options(ksp);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::attach_preconditioner(KSP &ksp) const {
        Impl w(ksp, false);
        auto delegate_ptr =
            std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());

        if (delegate_ptr) {
            if (ksp_->has_shell_pc()) {
                m_utopia_warning_once(
                    "set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
            }
        } else if (this->get_preconditioner()) {
            auto shell_ptr = this->get_preconditioner().get();
            w.attach_shell_preconditioner(UtopiaPCApplyShell, shell_ptr, nullptr, nullptr);
        }
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::set_preconditioner(const std::shared_ptr<Preconditioner> &precond) {
        PreconditionedSolver::set_preconditioner(precond);

        auto delegate_ptr =
            std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());

        if (delegate_ptr) {
            if (ksp_->has_shell_pc()) {
                m_utopia_warning_once(
                    "set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
                ksp_->pc_type("jacobi");
            }

        } else if (this->get_preconditioner()) {
            auto shell_ptr = this->get_preconditioner().get();
            ksp_->attach_shell_preconditioner(UtopiaPCApplyShell, shell_ptr, nullptr, nullptr);
        }
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::set_monitor_options(KSP &ksp) const {
        PetscErrorCode ierr;
        UTOPIA_UNUSED(ierr);
        if (this->verbose()) {
            ierr = KSPMonitorSet(ksp, MyKSPMonitor, new KSPLog(), [](void **arg) -> PetscErrorCode {
                auto ksp_log = (KSPLog **)(arg);
                delete *ksp_log;
                *ksp_log = nullptr;
                return 0;
            });
            assert(ierr == 0);
        }

        if (must_compute_cond_number()) {
            ierr = KSPSetComputeSingularValues(ksp, PETSC_TRUE);
            assert(ierr == 0);
        }
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::handle_reset(const Matrix &op) {
        bool must_reset = true;
        // bool must_reset = false;

        // if(this->has_operator()) {
        //     auto s = size(*this->get_operator());
        //     auto other_s = size(op);
        //     if(s != other_s) {
        //     must_reset = true;
        //     }
        // }

        // if(op.implementation().communicator() != ksp_->communicator())
        // {
        //     must_reset = true;
        // }

        if (must_reset) {
            auto temp_ksp = utopia::make_unique<Impl>(op.comm().get());
            temp_ksp->copy_settings_from(*ksp_);
            ksp_ = std::move(temp_ksp);

            if (this->get_preconditioner()) {
                set_preconditioner(this->get_preconditioner());
            }
        }
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op,
                                                  const std::shared_ptr<const Matrix> &prec) {
        UTOPIA_TRACE_REGION_BEGIN("KSPSolver::update(op,prec)");

        handle_reset(*op);
        set_monitor_options(ksp_->implementation());

        PreconditionedSolver::update(op, prec);
        ksp_->update(*op, *prec);

        UTOPIA_TRACE_REGION_END("KSPSolver::update(op,prec)");
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op) {
        UTOPIA_TRACE_REGION_BEGIN("KSPSolver::update");

        handle_reset(*op);

        PreconditionedSolver::update(op);
        set_monitor_options(ksp_->implementation());

        bool skip_set_operators = false;
        if (this->get_preconditioner()) {
            auto mat_prec =
                std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());
            if (mat_prec) {
                ksp_->update(*op, *mat_prec->get_matrix());
                skip_set_operators = true;
            }
        }

        if (!skip_set_operators) {
            ksp_->update(*op);
        }

        UTOPIA_TRACE_REGION_END("KSPSolver::update");
    }

    template <typename Matrix, typename Vector>
    KSPSolver<Matrix, Vector, PETSC> &KSPSolver<Matrix, Vector, PETSC>::operator=(
        const KSPSolver<Matrix, Vector, PETSC> &other) {
        if (this == &other) return *this;
        PreconditionedSolver::operator=(other);

        ksp_ = utopia::make_unique<Impl>(other.ksp_->communicator());
        ksp_->copy_settings_from(*other.ksp_);
        return *this;
    }

    template <typename Matrix, typename Vector>
    KSPSolver<Matrix, Vector, PETSC> &KSPSolver<Matrix, Vector, PETSC>::operator=(
        KSPSolver<Matrix, Vector, PETSC> &&other) {
        if (this == &other) return *this;
        PreconditionedSolver::operator=(std::move(other));
        ksp_ = std::move(other.ksp_);
        return *this;
    }

    template <typename Matrix, typename Vector>
    KSPSolver<Matrix, Vector, PETSC>::KSPSolver(const KSPSolver<Matrix, Vector, PETSC> &other)
        : PreconditionedSolverInterface<Vector>(other),
          PreconditionedSolver(other),
          ksp_(utopia::make_unique<Impl>(other.ksp_->communicator())) {
        ksp_->copy_settings_from(*other.ksp_);
    }

    template <typename Matrix, typename Vector>
    KSPSolver<Matrix, Vector, PETSC>::KSPSolver(KSPSolver<Matrix, Vector, PETSC> &&other)
        : PreconditionedSolver(std::move(other)), ksp_(std::move(other.ksp_)) {}

    template <typename Matrix, typename Vector>
    KSPSolver<Matrix, Vector, PETSC> *KSPSolver<Matrix, Vector, PETSC>::clone() const {
        return new KSPSolver(*this);
    }

    template <typename Matrix, typename Vector>
    typename utopia::KSPSolver<Matrix, Vector, PETSC>::Impl &KSPSolver<Matrix, Vector, PETSC>::ksp() {
        return *ksp_;
    }

    template <typename Matrix, typename Vector>
    const typename utopia::KSPSolver<Matrix, Vector, PETSC>::Impl &KSPSolver<Matrix, Vector, PETSC>::ksp() const {
        return *ksp_;
    }

    template <typename Matrix, typename Vector>
    KSP &KSPSolver<Matrix, Vector, PETSC>::implementation() {
        return ksp().implementation();
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::describe(std::ostream &os) const {
        os << "KSPSolver: \n";
        os << "ksp_type:       " << ksp_type() << "\n";
        os << "pc_type:        " << pc_type() << "\n";
        os << "solver_package: " << solver_package() << "\n";
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::read(Input &is) {
        PreconditionedSolver::read(is);

        // TODO add other options
        auto new_ksp_type = this->ksp_type();
        auto new_pc_type = this->pc_type();
        auto new_solver_package = this->solver_package();

        is.get("ksp_type", new_ksp_type);
        is.get("pc_type", new_pc_type);
        is.get("solver_package", new_solver_package);

        this->ksp_type(new_ksp_type);
        this->pc_type(new_pc_type);
        this->solver_package(new_solver_package);
    }

    template <typename Matrix, typename Vector>
    void KSPSolver<Matrix, Vector, PETSC>::print_usage(std::ostream &os) const {
        PreconditionedSolver::print_usage(os);
        this->print_param_usage(os, "ksp_type", "string", "Type of KSP solver.", "bcgs");
        this->print_param_usage(os, "pc_type", "string", "Type of preconditoner.", "jacobi");
        this->print_param_usage(os, "solver_package", "string", "Type of solver package.", " ");
    }
}  // namespace utopia

#endif  // UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP
