#ifndef UTOPIA_FE_BRATU_HPP
#define UTOPIA_FE_BRATU_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_libmesh.hpp"

#include <iostream>

namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector>
    class Bratu : public ExtendedFunction<Matrix, Vector> {
    public:
        using Scalar   = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using IndexSet = typename utopia::Traits<Vector>::IndexSet;

        Bratu(FunctionSpace &V, const Scalar lambda = 1.5) : V_(V), lambda_(lambda)

        {
            if(lambda>=0.0 && lambda < 6.81)
                lambda_=lambda;
            else{
                std::cout<<"Bratu: lambda is not in the correct range. Please choose from the range: 0 < lambda < 6.8.";
                lambda_ = 1.0;
            }

            initialize();
        }

        bool value(const Vector &x, Scalar &energy) const override
        {
            IndexSet ghost_nodes;
            convert(V_.dof_map().get_send_list(), ghost_nodes);
            Vector x_ = ghosted(V_.dof_map().n_local_dofs(), V_.dof_map().n_dofs(), ghost_nodes);
            x_ = x;
            synchronize(x_);

            auto u  = trial(V_);
            auto uk = interpolate(x_, u);

            auto f = 0.5 * inner(grad(uk), grad(uk)) * dX - lambda_ * exp(uk) * dX;
            utopia::assemble(f, energy);
            return true;
        }

        bool gradient_no_rhs(const Vector &x, Vector &gradient) const override
        {
            IndexSet ghost_nodes;
            convert(V_.dof_map().get_send_list(), ghost_nodes);
            Vector x_ = ghosted(V_.dof_map().n_local_dofs(), V_.dof_map().n_dofs(), ghost_nodes);
            x_ = x;
            synchronize(x_);

            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x_, u);

            auto l_form = inner(grad(uk), grad(v)) * dX - lambda_ * inner(exp( uk), v) * dX;

            utopia::assemble(l_form, gradient);

            apply_zero_boundary_conditions(V_.dof_map(), gradient);
            return true;
        }

        bool hessian(const Vector &x, Matrix &hessian) const override
        {
            IndexSet ghost_nodes;
            convert(V_.dof_map().get_send_list(), ghost_nodes);
            Vector x_ = ghosted(V_.dof_map().n_local_dofs(), V_.dof_map().n_dofs(), ghost_nodes);
            x_ = x;
            synchronize(x_);

            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x_, u);

            auto b_form = inner(grad(u), grad(v)) * dX - lambda_ * inner(exp(uk) * u, v) * dX;
            utopia::assemble(b_form, hessian);
            set_identity_at_constraint_rows(V_.dof_map(), hessian);
            return true;
        }

    private:
        FunctionSpace &V_;
        Scalar lambda_;

        void initialize()
        {
            Vector x = local_zeros(V_.dof_map().n_local_dofs());

            apply_boundary_conditions(V_.dof_map(), x);

            Vector marked;
            mark_constrained_dofs(V_.dof_map(), marked);
            this->set_equality_constrains(marked, x);
        }

    };

}

#endif // UTOPIA_FE_BRATU_HPP
