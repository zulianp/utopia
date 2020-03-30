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
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        typedef utopia::Preconditioner<Vector> Preconditioner;

        using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;

        using Super::solve;
        using Super::update;
        using Super::apply;

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

        void init_memory(const Layout &layout) override;

    private:
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
        Layout layout_;

    };
}

#endif //UTOPIA_BICG_STAB_HPP
