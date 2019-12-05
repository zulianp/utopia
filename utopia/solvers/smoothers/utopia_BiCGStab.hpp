#ifndef UTOPIA_BICG_STAB_HPP
#define UTOPIA_BICG_STAB_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Size.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_PreconditionedSolver.hpp"

namespace utopia {

    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class BiCGStab final : public OperatorBasedLinearSolver<Matrix, Vector>{
    public:
        typedef UTOPIA_SCALAR(Vector) 	 Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;

        using PreconditionedSolver<Matrix, Vector>::solve;

        BiCGStab();
        BiCGStab * clone() const override;

        inline bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override
        {
            if(this->get_preconditioner()) {
                return solve_preconditioned(A, b, x);
            } else {
                return solve_unpreconditioned(A, b, x);
            }
        }

        void update(const Operator<Vector> &A) override;


        void read(Input &in) override
        {
            OperatorBasedLinearSolver<Matrix, Vector>::read(in);
        }

        void print_usage(std::ostream &os) const override
        {
            OperatorBasedLinearSolver<Matrix, Vector>::print_usage(os);
        }

    private:
        void init(const SizeType &ls);
        bool solve_preconditioned(const Operator<Vector> &A, const Vector &b, Vector &x);
        bool solve_unpreconditioned(const Operator<Vector> &A, const Vector &b, Vector &x);

        Vector r0_;
        Vector r_;
        Vector v_;
        Vector p_;
        Vector h_;
        Vector y_;
        Vector t_;
        Vector s_;
        Vector z_;
        Vector K_inv_t_;

        bool initialized_;
        SizeType loc_size_;

    };
}

#endif //UTOPIA_BICG_STAB_HPP
