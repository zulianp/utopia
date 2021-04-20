#ifndef UTOPIA_PETSC_KSP_HPP
#define UTOPIA_PETSC_KSP_HPP

#include "utopia_Core.hpp"
#include "utopia_Input.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_petsc_Vector.hpp"

#include <petscksp.h>
#include <petscpc.h>
#include <petscsys.h>
#include <algorithm>

namespace utopia {

    // FIXME this whole wrapper needs fixing
    template <class Matrix, class Vector>
    class KSPWrapper;

    PetscErrorCode UtopiaPCApplyShell(PC pc, Vec x, Vec y);
    PetscErrorCode MyKSPMonitor(KSP, PetscInt, PetscReal, void *);
    std::string converged_str(KSPConvergedReason reason);

    /**@ingroup     Linear
     * @brief       Class provides interface to Petsc KSP solvers \n
     *              For setting up basic parameters, one can use classic Petsc runtime options, e.g.
     *              To see all possibilities, please refer to:
     *                                                   * <a
     * href="http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html">preconditioner types</a>
     *                                                   * <a
     * href="http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType">solver types</a>
     *
     *              Setting own/utopia preconditioner, can be done as following:
     *              \snippet tests/utopia_SolverTest.cpp PetscKSPSolver solve example1
     *              Detailed information about preconditioners, can be found in  \ref precondotioners.
     */
    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class KSPSolver {};

    template <typename Matrix, typename Vector>
    class KSPSolver<Matrix, Vector, PETSC> : public PreconditionedSolver<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Preconditioner = utopia::Preconditioner<Vector>;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "only works with petsc types");

        class Impl;

        KSPSolver();

        KSPSolver(std::unique_ptr<Impl> &&w);

        void wrap(KSP &ksp);

        ~KSPSolver() override;

        /* @brief      Sets the choice of direct solver.
         *             Please note, in petsc, direct solver is used as preconditioner alone, with proper settings.
         *
         * @param[in]  pc_type  The type of direct solver.
         */
        virtual void pc_type(const std::string &pc_type);

        /**
         * @brief      Sets KSP type
         */
        virtual void ksp_type(const std::string &ksp_type);

        /**
         * @brief      Sets solver package for choice of direct solver.
         *
         * @param[in]  package  The solver package.
         */
        void solver_package(const std::string &package);

        /**
         * @brief      Returns type of direct solver.
         */
        std::string pc_type() const;

        /**
         * @brief      Returns ksp package type
         */
        std::string ksp_type() const;

        /**
         * @brief      Setter for number of global subdomains
         */
        void number_of_subdomains(const SizeType &n);

        /**
         * @brief      Set norm type
         */
        void norm_type(const std::string &norm_type);

        /**
         * @brief      Setter for overlap used inside of Additive Schwarz method
         */
        void overlap(const SizeType &n);

        /**
         * @brief      Sets ksp and pc type for all sub_ksp and all sub_pc
         */
        void sub_ksp_pc_type(const std::string ksp_type, const std::string pc_type);

        /**
         * @brief      Sets ksp type for all sub_ksp
         */
        void sub_ksp_type(const std::string type);

        /**
         * @brief      Sets pc type for all sub_ksp
         */
        void sub_pc_type(const std::string type);

        /**
         * @brief      Sets solver package type for all sub_ksp
         */
        void sub_solver_package(const std::string type);

        /**
         * @brief      Returns type of solver package.
         */
        std::string solver_package() const;

        /**
         * @brief        Solve function for all PETS KSP solvers.
         *              It is also compatible with our own Utopia preconditioners.
         */
        bool apply(const Vector &b, Vector &x) override;

        bool smooth(const Vector &rhs, Vector &x) override;

        bool must_compute_cond_number() const;

        void describe(std::ostream &os) const;

        /**
         * @brief      Sets the default options for PETSC KSP solver. \n
         *             Default: BiCGstab
         *
         * @param      ksp   The ksp
         */
        virtual void set_ksp_options(KSP &ksp);

        virtual void attach_preconditioner(KSP &ksp) const;

        void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) override;

        virtual void set_monitor_options(KSP &ksp) const;

        void handle_reset(const Matrix &op);

        void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) override;

        void update(const std::shared_ptr<const Matrix> &op) override;

        KSPSolver &operator=(const KSPSolver &other);

        KSPSolver &operator=(KSPSolver &&other);

        KSPSolver(const KSPSolver &other);

        KSPSolver(KSPSolver &&other);

        KSPSolver *clone() const override;

        Impl &ksp();

        const Impl &ksp() const;

        void set_initial_guess_non_zero(const bool val);

        KSP &implementation();

        void read(Input &is) override;
        void print_usage(std::ostream &os = std::cout) const override;

        void factor_set_pivot_in_blocks(const bool val);

    protected:
        std::unique_ptr<Impl> ksp_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_KSP_HPP
