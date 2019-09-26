#ifndef UTOPIA_UI_FORCING_FUNCTION_HPP
#define UTOPIA_UI_FORCING_FUNCTION_HPP

#include "utopia_ForcedMaterial.hpp"
#include "utopia_ui.hpp"
#include "utopia_AdaptiveFunction.hpp"

namespace utopia {



    inline double surface_area(LibMeshFunctionSpace &V, int side)
    {
        auto &dof_map = V.dof_map();
        auto u = trial(V);
        UVector x = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
        x.set(1.);

        double area = -1.;

        utopia::assemble(
            surface_integral(interpolate(x, u), side),
            area
        );

        return area;
    }

    template<class FunctionSpace, class Vector>
    class UIForcingFunction final : public CompositeForcingFunction<Vector>, public Configurable {
    public:

        UIForcingFunction(FunctionSpace &V) : V_(V) {}

        ~UIForcingFunction() {}

        void read(Input &is) override {
            is.get_all([this](Input &is) {
                std::string clazz = "";

                is.get("class", clazz);

                if(clazz == "AdaptiveFunction") {
                    auto ff = std::const_pointer_cast<AdaptiveFunction<FunctionSpace,Vector>>(make_ref<AdaptiveFunction<FunctionSpace,Vector>>(V_));
                    //auto ff = make_ref<AdaptiveFunction<FunctionSpace,Vector>>(V_);
                    ff->read(is);
                    //ff->read(is);
                    this->add(ff);
                    std::cout << "AdaptiveFunction" << std::endl;
                    return;
                }


                int block = -1;

                std::string type = "volume";
                is.get("block", block);
                is.get("type", type);


                if(block == -1) {
                    std::cerr << "[Error]ForcingFunction block not specified" << std::endl;
                    return;
                }


#ifdef WITH_TINY_EXPR
                std::string value;
                is.get("value", value);
                auto f = symbolic(value);
#else
                double value = 0.;
                is.get("value", value);
                auto f = coeff(value);
#endif //WITH_TINY_EXPR

                if(type == "surface") {
                    int normalize_by_area = 0;
                    is.get("normalize-by-area", normalize_by_area);

                    double area = 1.;
                    if(normalize_by_area) {
                        area = surface_area(V_, block);
                        std::cout << "normalizing by area: " << area << std::endl;
                    }

                    auto v = test(V_);
                    auto l_form = surface_integral((1./area) * inner(f, v), block);

                    auto ff = std::make_shared<ConstantForcingFunction<Vector>>();
                    ff->init(l_form);
                    this->add(ff);
                } else {
                    auto v = test(V_);
                    auto l_form = integral(inner(f, v), block);

                    auto ff = std::make_shared<ConstantForcingFunction<Vector>>();
                    ff->init(l_form);
                    this->add(ff);
                }
            });
        }

    private:
        FunctionSpace &V_;
    };


    template<class FunctionSpace, class Vector>
    class UIForcingFunction<ProductFunctionSpace<FunctionSpace>, Vector> final : public CompositeForcingFunction<Vector>, public Configurable {
    public:

        UIForcingFunction(ProductFunctionSpace<FunctionSpace> &V) : V_(V) {}

        virtual ~UIForcingFunction() {}

        void read(Input &is) override {
            is.get_all([this](Input &is) {

                int block = -1;
                int coord = 0;

                std::string type = "volume";
                is.get("block", block);
                is.get("coord", coord);
                is.get("type", type);

#ifdef WITH_TINY_EXPR
                std::string value;
                is.get("value", value);
                auto f = symbolic(value);
#else
                double value = 0.;
                is.get("value", value);
                auto f = coeff(value);
#endif //WITH_TINY_EXPR

                if(type == "surface") {

                    int normalize_by_area = 0;
                    is.get("normalize-by-area", normalize_by_area);

                    double area = 1.;
                    if(normalize_by_area) {
                        area = surface_area(V_[0], block);
                        std::cout << "normalizing by area: " << area << std::endl;
                    }

                    auto v = test(V_[coord]);
                    auto l_form = surface_integral((1./area) * inner(f, v), block);

                    auto ff = std::make_shared<ConstantForcingFunction<Vector>>();
                    ff->init(l_form);
                    this->add(ff);
                } else {
                    auto v = test(V_[coord]);
                    auto l_form = integral(inner(f, v), block);

                    auto ff = std::make_shared<ConstantForcingFunction<Vector>>();
                    ff->init(l_form);
                    this->add(ff);
                }
            });
        }

    private:
        ProductFunctionSpace<FunctionSpace> &V_;
    };
}


#endif //UTOPIA_UI_FORCING_FUNCTION_HPP
