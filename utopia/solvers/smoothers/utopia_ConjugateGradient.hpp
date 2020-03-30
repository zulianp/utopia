#ifndef UTOPIA_CONJUGATE_GRADIENT_HPP
#define UTOPIA_CONJUGATE_GRADIENT_HPP

#include "utopia_MatrixFreeLinearSolver.hpp"
#include <memory>

namespace utopia {

    /**
     * @brief      Conjugate Gradient solver. Works with all utopia tensor types.
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class ConjugateGradient final : public OperatorBasedLinearSolver<Matrix, Vector>
    {
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::Preconditioner<Vector> Preconditioner;

    public:
        using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;
        using Super::solve;
        using Super::update;

        ConjugateGradient();

        void reset_initial_guess(const bool val);
        inline void apply_gradient_descent_step(const bool val)
        {
            apply_gradient_descent_step_ = val;
        }

        void read(Input &in) override;

        void init_memory(const Layout &layout) override;

        void print_usage(std::ostream &os) const override;

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override;

        void update(const Operator<Vector> &A) override;

        ConjugateGradient * clone() const override;

        void copy(const ConjugateGradient &other);
        ConjugateGradient(const ConjugateGradient &other);

    private:
        bool unpreconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x);

        bool preconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x);

        bool check_solution(const Operator<Vector> &A, const Vector &x, const Vector &b) const;

        void gradient_descent_step(
                const Operator<Vector> &A,
                const Vector &b,
                Vector &x
        );

        bool reset_initial_guess_;
        bool initialized_;
        bool apply_gradient_descent_step_;
        Layout layout_;

        //This fields are not to be copied anywhere
        Vector r, p, q, Ap, r_new, z, z_new;
    };
}


#endif //UTOPIA_CONJUGATE_GRADIENT_HPP

