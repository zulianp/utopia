#ifndef UTOPIA_ADAPTIVE_FUNCTION_HPP
#define UTOPIA_ADAPTIVE_FUNCTION_HPP

#include "utopia_ElasticMaterial.hpp"
#include "utopia_ForcedMaterial.hpp"
#include "utopia_fe_EDSL.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

namespace utopia {
    template <class FunctionSpace, class Vector>
    class AdaptiveFunction : public ForcingFunction<Vector>, public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        AdaptiveFunction(FunctionSpace &V) : V_(V), value_("0"), block_(-1), area_(1), type_("volume") {}

        ~AdaptiveFunction() {}

        void read(Input &is) override {
            // std::string type = "volume";
            is.get("block", block_);
            is.get("type", type_);

            if (block_ == -1) {
                std::cerr << "[Error] ForcingFunction block not specified" << std::endl;
                return;
            }

            is.get("value", value_);

            if (type_ == "surface") {
                int normalize_by_area = 0;
                is.get("normalize-by-area", normalize_by_area);

                double area = 1.;
                if (normalize_by_area) {
                    area_ = surface_area(V_, block_);
                    std::cout << "normalizing by area: " << area_ << std::endl;
                }
            }
        }

        bool eval(const Vector &, Vector &result) override {
#ifdef UTOPIA_ENABLE_TINY_EXPR
            auto f = symbolic(value_);
#else
            double value = atof(value_.c_str());
            auto f = coeff(value);
#endif  // UTOPIA_ENABLE_TINY_EXPR

            // std::cout<<"Ciao"<<std::endl;
            if (type_ == "surface") {
                auto v = test(V_);
                auto l_form = surface_integral((1. / area_) * inner(f, v), block_);
                assemble(l_form, result);
            } else {
                auto v = test(V_);
                auto l_form = integral(inner(f, v), block_);
                assemble(l_form, result);
            }

            // disp(result);

            return true;
        }

    private:
        FunctionSpace &V_;
        std::string value_;
        int block_;
        Scalar area_;
        std::string type_;
    };

}  // namespace utopia
#endif
