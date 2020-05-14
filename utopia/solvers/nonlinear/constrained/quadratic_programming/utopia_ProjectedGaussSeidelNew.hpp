#ifndef UTOPIA_PROJECTED_GAUSS_SEIDEL_NEW_HPP
#define UTOPIA_PROJECTED_GAUSS_SEIDEL_NEW_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_SolverForwardDeclarations.hpp"

#include <cmath>
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_ProjectedGaussSeidel.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {

    template <class Matrix, class LocalMat>
    inline void local_block_view(const Matrix &, LocalMat &) {
        static_assert(Traits<LocalMat>::Backend != PETSC, "IMPLEMENT ME");
        static_assert(Traits<LocalMat>::Backend == HOMEMADE, "IMPLEMENT ME");
    }

    // FIXME make it work for all other backends
    template <class Matrix, class Vector>
    class ProjectedGaussSeidel<Matrix, Vector, PETSC> : public QPSolver<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using Super = utopia::QPSolver<Matrix, Vector>;

        ~ProjectedGaussSeidel() override;
        ProjectedGaussSeidel();
        ProjectedGaussSeidel(const ProjectedGaussSeidel &other);

        ProjectedGaussSeidel *clone() const override;
        void read(Input &in) override;
        void print_usage(std::ostream &os) const override;
        bool smooth(const Vector &b, Vector &x) override;
        bool apply(const Vector &b, Vector &x) override;
        void init_memory(const Layout &layout) override;
        void update(const std::shared_ptr<const Matrix> &op) override;

        inline void use_line_search(const bool val) { use_line_search_ = val; }

        inline SizeType n_local_sweeps() const { return n_local_sweeps_; }

        inline void n_local_sweeps(const SizeType n_local_sweeps) { n_local_sweeps_ = n_local_sweeps; }

        inline void use_symmetric_sweep(const bool use_symmetric_sweep) { use_symmetric_sweep_ = use_symmetric_sweep; }
        inline void use_sweeper(const bool val) { use_sweeper_ = val; }

        inline void l1(const bool val) { l1_ = val; }

    protected:
        virtual bool step(const Matrix &A, const Vector &b, Vector &x);
        virtual bool unconstrained_step(const Matrix &A, const Vector &b, Vector &x);

    private:
        bool use_line_search_;
        bool use_symmetric_sweep_;
        bool l1_;

        SizeType n_local_sweeps_;
        SizeType check_s_norm_each_;

        Vector r, d, ub_loc, lb_loc, c, d_inv, x_old, descent_dir, Ac;
        Vector inactive_set_;
        Vector is_c_;

        bool use_sweeper_;
        std::unique_ptr<ProjectedGaussSeidelSweep<Scalar, SizeType> > sweeper_;

        void init(const Matrix &A);

        /// residual must be computed outside
        void apply_local_sweeps(const Matrix &A, const Vector &r, const Vector &lb, const Vector &ub, Vector &c) const;

        void non_linear_jacobi_step(const Matrix &A, const Vector &b, Vector &x);

        /// residual must be computed outside
        void apply_local_sweeps_unconstrained(const Matrix &A, const Vector &r, Vector &c) const;
    };
}  // namespace utopia

#endif  // UTOPIA_PROJECTED_GAUSS_SEIDEL_NEW_HPP
