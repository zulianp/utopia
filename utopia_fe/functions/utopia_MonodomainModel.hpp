#ifndef UTOPIA_MONODOMAIN_MODEL_HPP
#define UTOPIA_MONODOMAIN_MODEL_HPP

#include "utopia.hpp"
#include "utopia_fe_core.hpp"

namespace utopia {
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class MonodomainModel : public Function <Matrix, Vector, Backend> {
    public:
        DEF_UTOPIA_SCALAR(Matrix); 

        virtual ~MonodomainModel() { }

        virtual bool value(const Vector &point, Scalar &value) const
        {
            //FIXME
            return 0;
        }

        virtual bool gradient(const Vector &point, Vector &result) const
        {

        }

        virtual bool hessian(const Vector &x, Matrix &H) const
        {
            auto u = trial(*V_);
            auto v = test(*V_);

            if(empty(lapl_)) {
                auto b_form = (inner(grad(u), grad(v)) * dX);
                utopia::assemble(b_form, lapl_);
                lapl_ *= diffusivity_;
                // JacIonCurrent_ = lapl_;
            }

            // JacIonCurrent_ *= 0.;


            auto uk = interpolate(x, u);

            auto c1 = coeff(u_unst_ - u_max_ - u_min_);
            auto c2 = coeff(u_min_ * u_max_ + u_min_ * u_unst_ - u_max_ * u_unst_);
            auto j_form = alpha * inner(u * u - u * c1 + c2, v) * dX;
            utopia::assemble(j_form, JacIonCurrent_);

            H = lapl_ + JacIonCurrent_;
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
    Matrix JacIonCurrent_;
    Scalar diffusivity_;
    Scalar u_min_, u_max_, u_unst_;
    Scalar alpha_;
};
}
#endif //UTOPIA_MONODOMAIN_MODEL_HPP
