#ifndef UTOPIA_FLOW_POST_PROCESSOR_HPP
#define UTOPIA_FLOW_POST_PROCESSOR_HPP

#include "utopia_Traits.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_Coefficient.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include <memory>
#include <vector>
#include <fstream>

namespace utopia {

    template<class FunctionSpace, class Vector>
    class PostProcessor : public Configurable {
    public:
        virtual ~PostProcessor() {}
        virtual void apply(FunctionSpace &V, const Vector &sol) = 0;
        virtual void apply(FunctionSpace &V, const Vector &sol,  const Vector &other) = 0;
        virtual void describe(std::ostream &os = std::cout) const = 0;
        virtual void export_values() const = 0;
    };


    template<class FunctionSpace, class Vector>
    class FluxPostProcessor final : public PostProcessor<FunctionSpace, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        void read(Input &in) override 
        {
            in.get("side", side_);
            in.get("diffusivity", diffusivity_);

            filename_ += std::to_string(side_) + ".txt";
            in.get("file-name", filename_);

            assert(valid());
        }

        void apply(FunctionSpace &V, const Vector &pressure) override
        {   
            auto u = trial(V);
            auto v = test(V);

            auto p = interpolate(pressure, u);
            auto sampler_fun = ctx_fun(sampler_);

            auto vel = (sampler_fun * diffusion_tensor_) * grad(p);

            auto form_1 = surface_integral(
                inner( 
                    inner(vel, normal() ), v)
                , side_
            );

            auto form_2 = surface_integral(inner(coeff(1.0), v), side_);

            assemble(form_1, weighted_flux_);
            assemble(form_2, mass_vector_);

            weighted_flux_ *= -1.0;

            e_pseudo_inv(mass_vector_, mass_vector_, 1e-12);
            flux_ = e_mul(mass_vector_, weighted_flux_);

            write("flux_" + std::to_string(side_) + ".e", V, flux_);

            const Scalar min_flux = min(flux_);
            const Scalar max_flux = max(flux_);
            const Scalar int_flux = sum(weighted_flux_);

            std::cout << "flux(" << side_ << ") in [" << min_flux << ", " << max_flux << "], ";
            std::cout << "int(flux): " << int_flux << std::endl;
        }

        void apply(FunctionSpace &V, const Vector &pressure, const Vector &concentration) override
        {   
            // auto u = trial(V);
            // auto v = test(V);

            // auto p = interpolate(pressure, u);
            // auto c = interpolate(concentration, u);

            // auto sampler_fun = ctx_fun(sampler_);

            // auto vel = (sampler_fun * diffusion_tensor_) * grad(p);

            // auto form = surface_integral(
            //     inner( 
            //         inner(vel, normal() ), c)
            //     , side_
            // );

            // outflux_ = 0.0;
            // assemble(form, outflux_);
            // outflux_ = -outflux_;

            if(empty(flow_matrix_)) {

                auto c = trial(V);
                auto q = test(V);

                // auto sampler_fun = ctx_fun(sampler_);

                auto p = interpolate(pressure, c);
                auto vel = coeff(diffusivity_) * grad(p);
                auto flow_form = surface_integral(inner(c * inner(vel, normal()), q), side_);

                assemble(flow_form, flow_matrix_);
                flow_matrix_ *= -1.0;
            }

            outflux_ = sum(flow_matrix_ * concentration);

            all_outflux_.push_back(outflux_);
        }

        inline void sampler(const std::shared_ptr<UIFunction<double>> &s)
        {
            sampler_ = s;
        }

        inline void diffusion_tensor(const ElementMatrix &d)
        {
            diffusion_tensor_ = d;
        }

        inline Vector &flux() 
        {
            return flux_;
        }

        inline const Vector &flux() const
        {
            return flux_;
        }

        inline void side(const int side)
        {
            side_ = side;
        }

        bool valid() const
        {
            assert(side_ != -1);
            
            if(side_ == -1) {
                return false;
            }

            assert(bool(sampler_));

            if(!sampler_) return false;


            return true;
        }

        inline double outflux() const
        {
            return outflux_;
        }

        FluxPostProcessor()
        : side_(-1), outflux_(0.0), filename_("flux")
        {}

        void describe(std::ostream &os = std::cout) const override
        {
            os << "ouflux(" << side_ << ") = " << outflux_ << std::endl;
        }

        void export_values() const override
        {
            std::ofstream os(filename_.c_str());
            assert(os.good());

            if(os.good()) {
                for(auto f : all_outflux_) {
                    os << f << "\n";
                }
            }

            os.close();
        }

    private:
        Vector flux_;
        Vector weighted_flux_;
        Vector mass_vector_;
        double outflux_;
        double diffusivity_;
        std::shared_ptr<UIFunction<double>> sampler_;
        ElementMatrix diffusion_tensor_;
        int side_;
        USparseMatrix flow_matrix_;

        std::vector<double> all_outflux_;
        std::string filename_;
    };

}

#endif //UTOPIA_FLOW_POST_PROCESSOR_HPP
