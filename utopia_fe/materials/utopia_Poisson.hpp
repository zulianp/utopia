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
        
        Poisson(FunctionSpace &V) : V_(V)
        {
            initialize();
        }
        
        
        //HANDLE boundary conditions when computing energy
        bool value(const Vector &x, typename Vector::Scalar &energy) const override
        {
            // auto u  = trial(V_);
            // auto uk = interpolate(x, u);
            // auto r  = interpolate(rhs_, u);
            
            // utopia::assemble(0.5 * inner(grad(uk), grad(uk)) * dX, energy);
            // energy = dot(x, rhs_) - energy;

            energy = 0.5 * dot(H_ * x, x) - dot(x, rhs_);
            return true;
        }
        
        bool gradient(const Vector &x, Vector &gradient) const override
        {
            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x, u);
            
            auto l_form = inner(grad(uk), grad(v)) * dX;
            utopia::assemble(l_form, gradient);
            
            gradient -= rhs_;
            apply_zero_boundary_conditions(V_.dof_map(), gradient);
            return true;
        }
        
        bool hessian(const Vector &x, Matrix &hessian) const override
        {
            hessian = H_;
            return true;
        }
        
        bool update(const Vector &) override
        {
            return true;
        }
        
        
    private:
        FunctionSpace &V_;
        Vector rhs_;
        Matrix H_;
        
        void initialize()
        {
            
            auto u  = trial(V_);
            auto v  = test(V_);
            auto b_form = inner(grad(u), grad(v)) * dX;
            utopia::assemble(b_form, H_);
            set_identity_at_constraint_rows(V_.dof_map(), H_);
            
            utopia::assemble(inner(coeff(1.), v) * dX, rhs_);
            apply_boundary_conditions(V_.dof_map(), rhs_);
        }
        
    };
    
}

#endif // UTOPIA_POISSON_HPP
