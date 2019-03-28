#ifndef UTOPIA_MONODOMAIN_MODEL_HPP
#define UTOPIA_MONODOMAIN_MODEL_HPP

#include "utopia_libmesh_Types.hpp"
#include "utopia.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class MonodomainModel : public Function<Matrix, Vector> {
    public:
        DEF_UTOPIA_SCALAR(Matrix);

        virtual ~MonodomainModel() { }

        virtual bool value(const Vector &x, Scalar &value) const
        {
            //FIXME
            assert(false && "implement me");
            return -1;
        }

        virtual bool gradient(const Vector &x, Vector &result) const
        {
            auto u = trial(*V_);
            auto v = test(*V_);

            auto uk     = interpolate(x, u);
            auto u_min  = coeff(u_min_);
            auto u_max  = coeff(u_max_);
            auto u_unst = coeff(u_unst_);

            auto l = (uk - u_min) * (uk - u_max) * (uk - u_unst);
            auto f_form = alpha_ * inner(l, v) * dX;

            utopia::assemble(f_form, result);
            result -= lapl_ * x;
            return false;
        }

        virtual bool hessian(const Vector &x, Matrix &H) const
        {
            auto u = trial(*V_);
            auto v = test(*V_);

            auto uk     = interpolate(x, u);
            auto u_min  = coeff(u_min_);
            auto u_max  = coeff(u_max_);
            auto u_unst = coeff(u_unst_);

            auto l = ( (uk - u_min) * (uk - u_max) + (uk - u_min) * (uk - u_unst) + (uk - u_max) * (uk - u_unst) );
            auto j_form = alpha_ * inner(l * u, v) * dX;

            Matrix deriv_Ion;
            utopia::assemble(j_form, deriv_Ion);

            H = lapl_ - deriv_Ion;
            return true;
        }

        void initialize(const std::shared_ptr<FunctionSpace> &V)
        {
            V_ = V;

            auto u = trial(*V_);
            auto v = test(*V_);

            auto b_form = (inner(grad(u), grad(v)) * dX);
            utopia::assemble(b_form, lapl_);
            lapl_ *= diffusivity_;
        }

        void set_function_space(const std::shared_ptr<FunctionSpace> &V)
        {
            V_ = V;
        }

        void set_diffusivity(const Scalar diffusivity)
        {
            diffusivity_ = diffusivity;
        }

        void set_ion_current_params(const Scalar alpha, const Scalar u_min, const Scalar u_max, const Scalar u_unst)
        {
            alpha_ = alpha;
            u_min_ = u_min;
            u_max_ = u_max;
            u_unst_ = u_unst;
        }

    private:
        std::shared_ptr<FunctionSpace> V_;
        Matrix lapl_;
        Matrix deriv_Ion;
        Scalar diffusivity_;
        Scalar u_min_, u_max_, u_unst_;
        Scalar alpha_;
    };
}
#endif //UTOPIA_MONODOMAIN_MODEL_HPP
