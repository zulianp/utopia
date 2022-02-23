#ifndef UTOPIA_PETSC_PROJECTED_GAUSS_SEIDEL_HPP
#define UTOPIA_PETSC_PROJECTED_GAUSS_SEIDEL_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_SolverForwardDeclarations.hpp"

#include <cmath>
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_ProjectedGaussSeidel.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_Smoother.hpp"

#include "utopia_petsc_Types.hpp"

namespace utopia {

    template <>
    class ProjectedGaussSeidel<PetscMatrix, PetscVector, PETSC> : public QPSolver<PetscMatrix, PetscVector> {
    public:
        using Scalar = typename Traits<PetscVector>::Scalar;
        using SizeType = typename Traits<PetscVector>::SizeType;
        using Layout = typename Traits<PetscVector>::Layout;
        using Super = utopia::QPSolver<PetscMatrix, PetscVector>;

        using VariableBoundSolverInterface = Super::VariableBoundSolverInterface;
        using PreconditionedSolverInterface = Super::PreconditionedSolverInterface;

        ~ProjectedGaussSeidel() override;
        ProjectedGaussSeidel();
        ProjectedGaussSeidel(const ProjectedGaussSeidel &other);

        ProjectedGaussSeidel *clone() const override;
        void read(Input &in) override;
        void print_usage(std::ostream &os) const override;
        bool smooth(const PetscVector &b, PetscVector &x) override;
        bool apply(const PetscVector &b, PetscVector &x) override;
        void init_memory(const Layout &layout) override;
        void update(const std::shared_ptr<const PetscMatrix> &op) override;

        inline void use_line_search(const bool val) { use_line_search_ = val; }

        inline SizeType n_local_sweeps() const { return n_local_sweeps_; }

        inline void n_local_sweeps(const SizeType n_local_sweeps) { n_local_sweeps_ = n_local_sweeps; }

        inline void use_symmetric_sweep(const bool use_symmetric_sweep) { use_symmetric_sweep_ = use_symmetric_sweep; }
        inline void use_sweeper(const bool val) { use_sweeper_ = val; }

        inline void l1(const bool val) { l1_ = val; }

        void set_sweeper(std::unique_ptr<ProjectedGaussSeidelSweep<PetscMatrix> > &&sweeper);

    protected:
        virtual bool step(const PetscMatrix &A, const PetscVector &b, PetscVector &x);
        virtual bool unconstrained_step(const PetscMatrix &A, const PetscVector &b, PetscVector &x);

    private:
        bool use_line_search_{false};
        bool use_symmetric_sweep_{true};
        bool l1_{false};

        SizeType n_local_sweeps_;
        SizeType check_s_norm_each_;

        PetscVector r, d, ub_loc, lb_loc, c, d_inv, x_old, descent_dir, Ac;
        PetscVector inactive_set_;
        PetscVector is_c_;

        bool use_sweeper_{true};
        std::unique_ptr<ProjectedGaussSeidelSweep<PetscMatrix> > sweeper_;

        /// residual must be computed outside
        void apply_local_sweeps(const PetscMatrix &A,
                                const PetscVector &r,
                                const PetscVector &lb,
                                const PetscVector &ub,
                                PetscVector &c) const;

        void non_linear_jacobi_step(const PetscMatrix &A, const PetscVector &b, PetscVector &x);

        /// residual must be computed outside
        void apply_local_sweeps_unconstrained(const PetscMatrix &A, const PetscVector &r, PetscVector &c) const;

        UTOPIA_NVCC_PRIVATE
        void init(const PetscMatrix &A);
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_PROJECTED_GAUSS_SEIDEL_HPP