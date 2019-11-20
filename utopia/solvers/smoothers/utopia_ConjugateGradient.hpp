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
        typedef UTOPIA_SCALAR(Vector) 	 Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::Preconditioner<Vector> Preconditioner;

    public:
        using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;
        using Super::solve;
        using Super::update;

        ConjugateGradient();

        void reset_initial_guess(const bool val);

        void read(Input &in) override;

        void print_usage(std::ostream &os) const override;

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override;

        void update(const Operator<Vector> &A) override;

        ConjugateGradient * clone() const override;

    private:
        bool unpreconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x);

        bool preconditioned_solve(const Operator<Vector> &A, const Vector &b, Vector &x);

        bool check_solution(const Operator<Vector> &A, const Vector &x, const Vector &b) const;

        void init(const SizeType &ls);

        // std::shared_ptr<Preconditioner> precond_;
        Vector r, p, q, Ap, r_new, z, z_new;
        bool reset_initial_guess_;
        bool initialized_;
        SizeType loc_size_;
    };
}


#endif //UTOPIA_CONJUGATE_GRADIENT_HPP

