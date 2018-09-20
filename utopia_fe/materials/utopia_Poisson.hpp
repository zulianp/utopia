#ifndef UTOPIA_POISSON_HPP
#define UTOPIA_POISSON_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"

#include <vector>

namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector>
    class Poisson final : public Function<Matrix, Vector> {
    public:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;
        
        Poisson(FunctionSpace &V) : V_(V), rhs_value_(1.)
        {
            // initialize();
        }
        
        bool value(const Vector &x, typename Vector::Scalar &energy) const override
        {
            auto u  = trial(V_);
            auto uk = interpolate(x, u);
            auto f1 = 0.5 * inner(grad(uk), grad(uk)) * dX;
            auto f2 = inner(coeff(rhs_value_), uk) * dX;
            utopia::assemble(f2 - f1, energy);
            return true;
        }
        
        bool gradient(const Vector &x, Vector &gradient) const override
        {
            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x, u);

            // auto l_form = inner(coeff(rhs_value_), v) * dX -  inner(grad(uk), grad(v)) * dX;
            auto l_form = inner(grad(uk), grad(v)) * dX - inner(coeff(rhs_value_), v) * dX;
            utopia::assemble(l_form, gradient);
            apply_zero_boundary_conditions(V_.dof_map(), gradient);
            return true;
        }
        
        bool hessian(const Vector &x, Matrix &hessian) const override
        {
            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x, u);
            auto b_form = inner(grad(u), grad(v)) * dX;
            utopia::assemble(b_form, hessian);
            set_identity_at_constraint_rows(V_.dof_map(), hessian);
            return true;
        }

        bool update(const Vector &) override
        { 
            return true; 
        }
        
    private:
        FunctionSpace &V_;
        Scalar rhs_value_;

    };
    
}

#endif // UTOPIA_POISSON_HPP
