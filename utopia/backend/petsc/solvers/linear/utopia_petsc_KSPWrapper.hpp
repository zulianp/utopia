#ifndef UTOPIA_PETSC_KSP_WRAPPER_HPP
#define UTOPIA_PETSC_KSP_WRAPPER_HPP

#include "utopia_Preconditioner.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_SolutionStatus.hpp"

#include "utopia_Core.hpp"

#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>

namespace utopia {
    std::string converged_str(KSPConvergedReason reason);

    class KSPTypes final {
    public:
        inline static KSPTypes &instance()
        {
            static KSPTypes instance_;
            return instance_;
        }

        inline const std::string &ksp(const SizeType i)
        {
            return ksp_.at(i);
        }

        inline const std::string &pc(const SizeType i)
        {
            return pc_.at(i);
        }

        inline const std::string &package(const SizeType i)
        {
            return package_.at(i);
        }

        inline bool is_ksp_valid(const std::string &type) const
        {
            return in_array(type, ksp_);
        }

        inline bool is_pc_valid(const std::string &type) const
        {
            return in_array(type, pc_);
        }

        inline bool is_solver_package_valid(const std::string &type) const
        {
            return in_array(type, package_);
        }

    private:
        KSPTypes();
        std::vector<std::string> ksp_;              /*!< Valid options for direct solver types. */
        std::vector<std::string> pc_;
        std::vector<std::string> package_;
    };

    template<class Matrix, class Vector>
    class KSPWrapper final {
    public:
        KSPWrapper(MPI_Comm comm,
                   const Parameters &params = Parameters())
        : ksp_(nullptr), owner_(true)
        {
            init(comm);

            ksp_type(KSPTypes::instance().ksp(0));
            pc_type(KSPTypes::instance().pc(0));
            solver_package(KSPTypes::instance().package(0));
        }

        KSPWrapper(KSP &ksp, const bool owner = false)
        : ksp_(ksp), owner_(owner)
        {}

        /* @brief      Sets the choice of direct solver.
         *             Please note, in petsc, direct solver is used as preconditioner alone, with proper settings.
         *
         * @param[in]  PCType  The type of direct solver.
         */
        inline void pc_type(const std::string &pc_type)
        {
            if(KSPTypes::instance().is_pc_valid(pc_type)) {
                PC pc;
                KSPGetPC(ksp_, &pc);
                PCSetType(pc, pc_type.c_str());
            }
        }

        /**
         * @brief      Sets KSP type
         */
        inline void ksp_type(const std::string & ksp_type)
        {
            if(KSPTypes::instance().is_ksp_valid(ksp_type)) {
                KSPSetType(ksp_, ksp_type.c_str());
            }
        }

        /**
         * @brief      Sets solver package for choice of direct solver.
         *
         * @param[in]  SolverPackage  The solver package.
         */
        inline void solver_package(const std::string &package)
        {
            if(KSPTypes::instance().is_solver_package_valid(package)) {
#if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
                PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
                PC pc;
                ierr = KSPGetPC(ksp_, &pc);  assert(ierr == 0);
                ierr = PCFactorSetMatSolverPackage(pc, package.c_str()); assert(ierr == 0);
#else
                m_utopia_warning_once("PCFactorSetMatSolverPackage not available in petsc 3.9.0 find equivalent?");
#endif
            }
        }


        /**
         * @brief      Returns type of direct solver.
         */
        inline PCType pc_type() const
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            PC pc;
            PCType ret;
            ierr = KSPGetPC(ksp_, &pc); assert(ierr == 0);
            ierr = PCGetType(pc, &ret); assert(ierr == 0);
            return ret;
        }

        inline bool has_shell_pc() const
        {
            return std::string(PCSHELL) == pc_type();
        }

        /**
         * @brief      Returns ksp package type
         */
        inline KSPType ksp_type() const
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            KSPType ret;
            ierr = KSPGetType(ksp_, &ret); assert(ierr == 0);
            return ret;
        }

        /**
         * @brief      Returns type of solver package.
         */
        inline std::string solver_package() const
        {
#if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            const MatSolverPackage stype;
            PC pc;
            ierr = KSPGetPC(ksp_, &pc);                     assert(ierr == 0);
            ierr = PCFactorGetMatSolverPackage(pc, &stype); assert(ierr == 0);
            return stype;
#else
            static const char * ret = " ";
            return ret;
#endif
        }

        inline  bool is_null() const
        {
            return ksp_ == nullptr;
        }

        ~KSPWrapper()
        {
            destroy();
        }

        inline void destroy()
        {
            if(ksp_) {
                if(owner_) {
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

        void describe() const
        {
            PetscBool bool_value;
            KSPGetInitialGuessNonzero(ksp_, &bool_value);
            std::cout << "-------------------------------------------" << std::endl;
            std::cout << "KSPGetInitialGuessNonzero: " << bool_value << std::endl;
            std::cout << "PCType:                    " << pc_type() << std::endl;
            std::cout << "KSPType:                   " << ksp_type() << std::endl;
            std::cout << "-------------------------------------------" << std::endl;

            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            ierr = KSPView(ksp_, PETSC_VIEWER_STDOUT_(communicator())); assert(ierr == 0);
        }

        inline void init(const MPI_Comm comm)
        {
            destroy();
            KSPCreate(comm, &ksp_);
            KSPSetFromOptions(ksp_);
            KSPSetComputeSingularValues(ksp_, PETSC_FALSE);
            // KSPSetNormType(ksp_, KSP_NORM_UNPRECONDITIONED);
        }

        inline void set_initial_guess_non_zero(const bool val)
        {
            if(val) {
                KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
            } else {
                KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);
            }
        }

        void pc_copy_settings(PC &other_pc) const
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            PC this_pc;
            PCType type;

            ierr = KSPGetPC(ksp_, &this_pc);    assert(ierr == 0);

            ierr = PCGetType(this_pc, &type);   assert(ierr == 0);
            ierr = PCSetType(other_pc, type);   assert(ierr == 0);

#if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
            const MatSolverPackage stype;
            ierr = PCFactorGetMatSolverPackage(this_pc, &stype); assert(ierr == 0);
            ierr = PCFactorSetMatSolverPackage(other_pc, stype); assert(ierr == 0);
#endif
        }

        void copy_settings_from(const KSPWrapper &other)
        {
            other.copy_settings_to(ksp_);
        }

        void copy_settings_to(KSP &other_ksp) const
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            PetscBool bool_value;
            KSPType type;
            PetscReal rtol;
            PetscReal abstol;
            PetscReal dtol;
            PetscInt maxits;
            PC other_pc;

            KSPSetFromOptions(other_ksp);

            ierr = KSPGetComputeSingularValues(ksp_, &bool_value);     assert(ierr == 0);
            ierr = KSPSetComputeSingularValues(other_ksp, bool_value); assert(ierr == 0);

            ierr = KSPGetType(ksp_, &type);                            assert(ierr == 0);
            ierr = KSPSetType(other_ksp, type);                        assert(ierr == 0);

            ierr = KSPGetInitialGuessNonzero(ksp_, &bool_value);       assert(ierr == 0);
            ierr = KSPSetInitialGuessNonzero(other_ksp, bool_value);   assert(ierr == 0);

            ierr = KSPGetTolerances(ksp_, &rtol, &abstol, &dtol, &maxits);  assert(ierr == 0);
            ierr = KSPSetTolerances(other_ksp, rtol, abstol, dtol, maxits); assert(ierr == 0);

            ierr = KSPGetPC(other_ksp, &other_pc);

            pc_copy_settings(other_pc);
        }

        void attach_shell_preconditioner(
            PetscErrorCode (*apply)(PC, Vec, Vec),
            void *ctx,
            PetscErrorCode (*setup)(PC),
            PetscErrorCode (*destroy)(PC)
            )
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            PC pc;

            ierr = KSPGetPC(ksp_, &pc);                         assert(ierr == 0);
            ierr = PCSetType(pc, PCSHELL);                      assert(ierr == 0);
            ierr = PCShellSetApply(pc, apply);                  assert(ierr == 0);

            ierr = PCShellSetName(pc, "Utopia Preconditioner"); assert(ierr == 0);

            if(setup) {
                ierr = PCShellSetSetUp(pc, setup);              assert(ierr == 0);
            }

            if(destroy) {
                ierr = PCShellSetSetUp(pc, destroy);            assert(ierr == 0);
            }

            if(ctx) {
                ierr = PCShellSetContext(pc, ctx);              assert(ierr == 0);
            }

            //usefull methods
            //PCGetOperatorsSet(PC pc,PetscBool  *mat,PetscBool  *pmat)
            //PetscErrorCode  PCGetOperators(PC pc,Mat *Amat,Mat *Pmat)
            //PetscErrorCode  PCShellSetApplyRichardson(PC pc,PetscErrorCode (*apply)(PC,Vec,Vec,Vec,PetscReal,PetscReal,PetscReal,PetscInt,PetscBool,PetscInt*,PCRichardsonConvergedReason*))
            //PCGetReusePreconditioner(PC pc,PetscBool *flag)
        }

        void update(const Matrix &mat)
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);

            ierr = KSPSetOperators(ksp_, raw_type(mat), raw_type(mat)); assert(ierr == 0);
            ierr = KSPSetUp(ksp_);                                      assert(ierr == 0);
        }

        void update(const Matrix &mat, const Matrix &prec)
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);

            ierr = KSPSetOperators(ksp_, raw_type(mat), raw_type(prec)); assert(ierr == 0);
            ierr = KSPSetUp(ksp_);                                       assert(ierr == 0);
        }

        bool smooth(const SizeType sweeps,
                    const Vector &rhs,
                    Vector &x)
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            KSPNormType normtype;
            PetscReal rtol;
            PetscReal abstol;
            PetscReal dtol;
            PetscInt maxits;
            void *ctx;

            //back-up settings
            ierr = KSPGetNormType(ksp_, &normtype);                                                    assert(ierr == 0);
            ierr = KSPGetTolerances(ksp_, &rtol, &abstol, &dtol, &maxits);                             assert(ierr == 0);

            //set up smoothing things
            ierr = KSPSetTolerances(ksp_, 0., 0., PETSC_DEFAULT, sweeps);                              assert(ierr == 0);
            ierr = KSPSetNormType(ksp_, KSP_NORM_NONE);                                                assert(ierr == 0);
            ierr = KSPSetConvergenceTest(ksp_, KSPConvergedSkip, NULL, NULL);                          assert(ierr == 0);

            //perform smoothing
            ierr = KSPSolve(ksp_, raw_type(rhs), raw_type(x));                                         assert(ierr == 0);

            //reset-solver-stuff
            ierr = KSPSetNormType(ksp_, normtype);                                                     assert(ierr == 0);
            ierr = KSPSetTolerances(ksp_, rtol, abstol, dtol, maxits);                                 assert(ierr == 0);
            ierr = KSPConvergedDefaultCreate(&ctx);                                                    assert(ierr == 0);
            ierr = KSPSetConvergenceTest(ksp_,  KSPConvergedDefault, ctx, KSPConvergedDefaultDestroy); assert(ierr == 0);
            return true;
        }

        void set_tolerances(const PetscReal rtol,
                            const PetscReal atol,
                            const PetscReal dtol,
                            const PetscInt max_it)
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            ierr = KSPSetTolerances(ksp_, rtol, atol, dtol, max_it); assert(ierr == 0);
        }

        void set_monitor(
            PetscErrorCode (*monitor)(KSP,PetscInt,PetscReal,void *),
            void *mctx,
            PetscErrorCode (*monitordestroy)(void**)
        )
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            ierr = KSPMonitorSet(ksp_, monitor, mctx, monitordestroy); assert(ierr == 0);
        }

        bool apply(const Vector &b, Vector &x)
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            auto ls = local_size(b).get(0);
            auto gs = size(b).get(0);

            if(gs != size(x).get(0))
            {
                x = local_zeros(ls);
            }

            assert(ls == local_size(x).get(0));

            KSPConvergedReason  reason;
            ierr = KSPSolve(ksp_, raw_type(b), raw_type(x)); assert(ierr == 0);
            ierr = KSPGetConvergedReason(ksp_, &reason);     assert(ierr == 0);

            if(reason < 0) {
                utopia_warning("ksp apply returned " + std::to_string(reason) + " = " + converged_str(reason) );
            }

            return reason >= 0;
        }

        void solution_status(SolutionStatus &status)
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
            PetscInt its;
            KSPConvergedReason reason;
            PetscReal rnorm;

            ierr = KSPGetIterationNumber(ksp_, &its);    assert(ierr == 0);
            ierr = KSPGetConvergedReason(ksp_, &reason); assert(ierr == 0);
            ierr = KSPGetResidualNorm(ksp_, &rnorm);     assert(ierr == 0);

            status.iterates = its;
            status.reason = reason;
            status.gradient_norm = rnorm;
        }

        void set_parameters(const Parameters params)
        {
            // our param options - useful for passing from passo
            ksp_type(params.lin_solver_type());
            pc_type(params.preconditioner_type());
            solver_package(params.preconditioner_factor_mat_solver_package());
            set_tolerances(params.rtol(), params.atol(), PETSC_DEFAULT, params.max_it());

            // petsc command line options
            char           name_[1024];
            PetscBool      flg;

#if UTOPIA_PETSC_VERSION_LESS_THAN(3,7,0)
            PetscOptionsGetString(nullptr, "-ksp_type", name_, 1024, &flg);
            if(flg) {
                ksp_type(name_);
            }

            PetscOptionsGetString(nullptr, "-pc_type", name_, 1024, &flg);
            if(flg) {
                pc_type(name_);
            }

            PetscOptionsGetString(nullptr, "-pc_factor_mat_solver_package", name_, 1024, &flg);
            if(flg) {
                solver_package(name_);
            }
#else
            PetscOptionsGetString(nullptr, nullptr, "-ksp_type", name_, 1024, &flg);
            if(flg) {
                ksp_type(name_);
            }

            PetscOptionsGetString(nullptr, nullptr, "-pc_type", name_, 1024, &flg);
            if(flg) {
                pc_type(name_);
            }

            PetscOptionsGetString(nullptr, nullptr, "-pc_factor_mat_solver_package", name_, 1024, &flg);
            if(flg) {
                solver_package(name_);
            }
#endif
        }

        KSP &implementation()
        {
            return ksp_;
        }

        const KSP &implementation() const
        {
            return ksp_;
        }

    private:
        KSP            ksp_;
        bool owner_;
    };

}



#endif //UTOPIA_PETSC_KSP_WRAPPER_HPP
