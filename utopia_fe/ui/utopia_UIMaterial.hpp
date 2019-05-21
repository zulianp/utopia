#ifndef UTOPIA_UI_MATERIAL_FUNCTION_HPP
#define UTOPIA_UI_MATERIAL_FUNCTION_HPP

#include "utopia_ElasticMaterial.hpp"
#include "utopia_StabilizedMaterial.hpp"
#include "utopia_LameeParameters.hpp"
#include "utopia_LinearElasticity.hpp"
#include "utopia_NeoHookean.hpp"
#include "utopia_SaintVenantKirchoff.hpp"

#include "utopia_ui.hpp"

namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector>
    class UIMaterial final : public ElasticMaterial<Matrix, Vector>, public Configurable {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        UIMaterial(FunctionSpace &V) : V_(V) {}

        ~UIMaterial() {}

        void read(Input &is) override {

            std::string material = "LinearElasticity";
            std::string stabilization = "none";
            Scalar stabilization_mag = 0.0001;

            is.get("material", material);
            is.get("stabilization", stabilization);
            is.get("stabilization-mag", stabilization_mag);
            is.get("parameters", params);

            params.describe(std::cout);

            if(material == "NeoHookean") {
                std::cout << "Using: NeoHookean" << std::endl;
                material_ = std::make_shared<NeoHookean<FunctionSpace, Matrix, Vector>>(V_, params);
            } else if(material == "SaintVenantKirchoff") {
                std::cout << "Using: SaintVenantKirchoff" << std::endl;
                material_ = std::make_shared<SaintVenantKirchoff<FunctionSpace, Matrix, Vector>>(V_, params);
            } else /*if(material == "LinearElasticity")*/ {
                material_ = std::make_shared<LinearElasticity<FunctionSpace, Matrix, Vector>>(V_, params);
            }

            if(stabilization != "none") {
                std::cout << "using stabilization: " << stabilization << " mag: " << stabilization_mag << std::endl;
                // StabilizedMaterial<FunctionSpace, Matrix, Vector> sm(V_, stabilization_mag, material_, stabilization);
                material_ = std::make_shared<StabilizedMaterial<FunctionSpace, Matrix, Vector>>(V_, stabilization_mag, material_, stabilization);
            }
        }

        inline bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            assert(material_);
            return material_->assemble_hessian_and_gradient(x, hessian, gradient);
        }

        inline bool stress(const Vector &x, Vector &result) override {
            assert(material_);
            return material_->stress(x, result);
        }

        inline void clear() override {
            assert(material_);
            material_->clear();
        }

        inline bool good() const
        {
            return static_cast<bool>(material_);
        }

        inline bool is_linear() const override { assert(material_); return material_->is_linear(); }

        LameeParameters params;

    private:
        FunctionSpace &V_;
        std::shared_ptr<ElasticMaterial<Matrix, Vector>> material_;
    };
}


#endif //UTOPIA_UI_MATERIAL_FUNCTION_HPP
