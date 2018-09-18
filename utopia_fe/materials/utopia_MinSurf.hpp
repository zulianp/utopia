#ifndef UTOPIA_MIN_SURF_HPP
#define UTOPIA_MIN_SURF_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"

#include <vector>

namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector>
    class MinSurf : public ExtendedFunction<Matrix, Vector> {
    public:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;
        
        MinSurf(FunctionSpace &V) : V_(V)
        {
            initialize();
        }
        
        bool value(const Vector &x, typename Vector::Scalar &energy) const override
        {
            // auto u  = trial(V_);
            // auto uk = interpolate(x, u);
            // auto f = sqrt( coeff(1) + inner(grad(uk), grad(uk)) * dX );
            // utopia::assemble(f, energy);
            return true;
        }
        
        bool gradient_no_rhs(const Vector &x, Vector &gradient) const override
        {
            // auto u  = trial(V_);
            // auto v  = test(V_);
            // auto uk = interpolate(x, u);
            // auto denom = sqrt(coeff(1) + inner(grad(uk), grad(uk)));
            // auto num   = inner(grad(uk), grad(v));
            // auto l_form = (num/denom) * dX;
            // utopia::assemble(l_form, gradient);
            return true;
        }
        
        bool hessian(const Vector &x, Matrix &hessian) const override
        {
            // auto u  = trial(V_);
            // auto v  = test(V_);
            // auto uk = interpolate(x, u);

            // auto denom = sqrt(coeff(1) + inner(grad(uk), grad(uk)));
            // auto b_form = grad(v)/denom -  inner(grad(v), inner(grad(uk), grad(u)) * grad(uk) * pow(denom, -3.) );
            // utopia::assemble(b_form, hessian);
            // set_identity_at_constraint_rows(V_.dof_map(), hessian);
            return true;
        }
        
    private:
        FunctionSpace &V_;

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

#endif // UTOPIA_MIN_SURF_HPP
