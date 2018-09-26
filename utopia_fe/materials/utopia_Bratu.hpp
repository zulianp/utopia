#ifndef UTOPIA_FE_BRATU_HPP
#define UTOPIA_FE_BRATU_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_libmesh.hpp"


namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector>
    class Bratu : public ExtendedFunction<Matrix, Vector> {
    public:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;
        
        Bratu(FunctionSpace &V, const Scalar lambda = 0.) : V_(V), lambda_(lambda)
        {
            initialize();
        }
        
        bool value(const Vector &x, Scalar &energy) const override
        {
            Vector x_ = ghosted(V_.dof_map().n_local_dofs(), V_.dof_map().n_dofs(), V_.dof_map().get_send_list()); 
            x_ = x;
            synchronize(x_);

            auto u  = trial(V_);
            auto uk = interpolate(x_, u);

            auto f = 0.5 * inner(grad(uk), grad(uk)) * dX - exp(lambda_ * uk) * dX;
            utopia::assemble(f, energy);
            return true;
        }
        
        bool gradient_no_rhs(const Vector &x, Vector &gradient) const override
        {
            Vector x_ = ghosted(V_.dof_map().n_local_dofs(), V_.dof_map().n_dofs(), V_.dof_map().get_send_list()); 
            x_ = x;
            synchronize(x_);

            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x_, u);

            auto l_form = inner(grad(uk), grad(v)) * dX - inner(exp(lambda_ * uk), v) * dX; 
            utopia::assemble(l_form, gradient);
            apply_zero_boundary_conditions(V_.dof_map(), gradient);
            return true;
        }
 
        bool hessian(const Vector &x, Matrix &hessian) const override
        {
            Vector x_ = ghosted(V_.dof_map().n_local_dofs(), V_.dof_map().n_dofs(), V_.dof_map().get_send_list()); 
            x_ = x;
            synchronize(x_);

            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x_, u);

            auto b_form = inner(grad(u), grad(v)) * dX - inner(exp(lambda_ * uk) * u, v) * dX;
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
