#ifndef UTOPIA_UI_MATERIAL_FUNCTION_HPP
#define UTOPIA_UI_MATERIAL_FUNCTION_HPP

#include "utopia_ElasticMaterial.hpp"
#include "utopia_LameeParameters.hpp"
#include "utopia_ui.hpp"

namespace utopia {
    template<class FunctionSpace, class Matrix, class Vector>
    class UIMaterial final : public ElasticMaterial<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        UIMaterial(FunctionSpace &V) : V_(V) {}

        ~UIMaterial() {}

        void read(Input &is) override;

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
