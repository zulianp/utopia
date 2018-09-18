#ifndef UTOPIA_BRATU_HPP
#define UTOPIA_BRATU_HPP

#include "utopia_Core.hpp"
#include <vector>

namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector>
    class Bratu : public ExtendedFunction<Matrix, Vector> {
    public:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;
        
        Bratu(FunctionSpace &V, const Scalar lambda = 0.2) : V_(V), lambda_(lambda)
        {
            initialize();
        }
        
        bool value(const Vector &x, typename Vector::Scalar &energy) const override
        {
            // auto u  = trial(V_);
            // auto uk = interpolate(x, u);

            // auto f = 0.5 * inner(grad(uk), grad(uk)) * dX - exp(lambda_ * uk) * dX;
            // utopia::assemble(f, energy);
            return true;
        }
        
        bool gradient_no_rhs(const Vector &x, Vector &gradient) const override
        {
            // auto u  = trial(V_);
            // auto v  = test(V_);
            // auto uk = interpolate(x, u);

            // auto l_form = inner(grad(uk), grad(v)) * dX - inner(exp(lambda_ * uk), v) * dX; 
            // utopia::assemble(l_form, gradient);
            return true;
        }
 
        bool hessian(const Vector &x, Matrix &hessian) const override
        {
            // auto u  = trial(V_);
            // auto v  = test(V_);
            // auto uk = interpolate(x, u);

            // auto b_form = inner(grad(u), grad(v)) * dX - inner(exp(lambda_ * uk) * u, v) * dX;
            // utopia::assemble(b_form, hessian);
            // set_identity_at_constraint_rows(V_.dof_map(), hessian);
            return true;
        }
        
    private:
        FunctionSpace &V_;
        Scalar lambda_;

        void initialize()
        {
            // Vector x = local_zeros(V_.dof_map().n_local_dofs());
            // apply_boundary_conditions(V_.dof_map(), x);

            // Vector marked;
            // mark_constrained_dofs(V_.dof_map(), marked);
            // this->set_equality_constrains(marked, x);
        }
        
    };
    
}

#endif // UTOPIA_BRATU_HPP
