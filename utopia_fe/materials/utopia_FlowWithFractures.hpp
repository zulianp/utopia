#ifndef UTOPIA_FLOW_WITH_FRACTURES_HPP
#define UTOPIA_FLOW_WITH_FRACTURES_HPP

#include "utopia_Model.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Integral.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_Flow.hpp"

#include <iostream>

namespace utopia {

    template<class FunctionSpace, class Matrix, class Vector>
    class FlowWithFractures final : public Model<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);
        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        typedef typename TraitsT::Vector ElementVector;

        FlowWithFractures(FunctionSpace &space)
        : space_(space), rescale_(1.0)
        {}

        class Fracture : public Configurable {
        public:
            Scalar angle;
            Scalar aperture;
            Scalar length;
            Scalar permeability;

            ElementVector center;
            ElementVector normal;
            ElementMatrix tangents;

            Fracture(const SizeType dim) 
            : angle(1.0), aperture(1.0), length(1.0), permeability(1.0)
            {
                center   = zeros(dim);
                normal   = zeros(dim);
                tangents = zeros(dim-1, dim);
            }

            static std::string coord_str(const SizeType d)
            {
                static const std::vector<std::string> str = {
                   "x", "y", "z", "t"
                };

                return str[d];
            }

            void read(Input &in) override
            {
                in.get("angle", angle);
                in.get("aperture", aperture);
                in.get("length", length);
                in.get("permeability", permeability);

                const SizeType n = size(center); 

                for(SizeType i = 0; i < n; ++i) {
                    Scalar v = 0.0;
                    in.get("center-" + coord_str(i), v);
                    center.set(i, v);

                    // v = 0.0;
                    // in.get("normal-" + coord_str(i), v);
                    // normal.set(i, v);

                    for(SizeType k = 0; k < n-1; ++k) {
                        v = 0.0;
                        in.get("tangent-" + std::to_string(k) + "-" + coord_str(i), v);
                        tangents.set(k, i, v);
                    }
                }

                if(n == 2) {
                    normal.set(0, -tangents.get(0, 1));
                    normal.set(1,  tangents.get(0, 0));
                } else 
                // if(n == 3) 
                {
                    //TODO
                    assert(false);
                    // ElementVector u = zeros(2), v = zeros(2);

                    // u.set(0, tangents.get(0, 1));
                    // u.set(1, tangents.get(0, 0));

                    // v.set(0, tangents.get(1, 1));
                    // v.set(1, tangents.get(1, 0));

                    // normal = cross(u, v);
                    // normal /= norm2(normal);
                }
            }
        };

        inline bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        {
            Chrono c;
            c.start();

            auto u = trial(space_);
            auto v = test(space_);

            auto bilinear_form = inner(
                    constant_permeability_ * grad(u), 
                    ctx_fun(permeability_.sampler()) * grad(v)
                ) * dX;

            utopia::assemble(bilinear_form, hessian);
            
            if(rescale_ != 1.0) {
                hessian *= rescale_;
                gradient *= rescale_;
            }

            c.stop();
            std::cout << "Flow assemly time: " << c << std::endl;

            const Scalar sum_g = sum(gradient);
            std::cout << "sum_g: " << sum_g << std::endl;
            return true;
        }

        inline bool is_linear() const override { return true; }

        inline void clear() override {}

        inline void read(Input &in) override 
        {
            int dim = space_.mesh().spatial_dimension();

            read_permeability_tensor(in);
            in.get("permeability-function", permeability_);
            in.get("fractures", [this, dim](Input &in) {
                in.get_all([this, dim](Input &in) {
                    auto f = utopia::make_unique<Fracture>(dim);
                });
            });
        }

        inline void rescale(const Scalar rescale)
        {
            rescale_ = rescale;
        }

    private:
        FunctionSpace &space_;
        UIScalarFunction<Scalar> permeability_;
        Scalar constant_permeability_;
        Scalar rescale_;
        std::vector<std::unique_ptr<Fracture>> fractures_;

        void read_permeability_tensor(Input &in)
        {
            constant_permeability_ = 1.;
            in.get("permeability",   constant_permeability_);
        }
    };

}

#endif //UTOPIA_FLOW_WITH_FRACTURES_HPP
