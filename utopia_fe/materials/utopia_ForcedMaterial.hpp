#ifndef UTOPIA_FORCED_MATERIAL_HPP
#define UTOPIA_FORCED_MATERIAL_HPP

#include "utopia.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

namespace utopia {

    template <class Vector>
    class ForcingFunction {
    public:
        virtual ~ForcingFunction() {}

        virtual bool eval(const Vector &x, Vector &result) {
            result = x;
            result.set(0.);
            return true;
        }
    };

    template <class Vector>
    class ConstantForcingFunction : public ForcingFunction<Vector> {
    public:
        bool eval(const Vector &, Vector &result) override {
            result = value_;
            return true;
        }

        template <class LinearForm>
        void init(const LinearForm &linear_form) {
            utopia::assemble(linear_form, value_);
        }

        void set(const Vector &value) { this->value_ = value; }

        Vector &value() { return value_; }

        inline const Vector &value() const { return value_; }

    private:
        Vector value_;
    };

    template <class Vector>
    class CompositeForcingFunction : public ForcingFunction<Vector> {
    public:
        virtual ~CompositeForcingFunction() {}

        bool eval(const Vector &x, Vector &result) override {
            result = x;
            result.set(0.);

            for (auto f_ptr : functions_) {
                if (!empty(buff_)) {
                    buff_.set(0.);
                }

                if (!f_ptr->eval(x, buff_)) {
                    assert(false);
                    return false;
                }

                result += buff_;
            }

            return true;
        }

        void add(const std::shared_ptr<ForcingFunction<Vector>> &function) { functions_.push_back(function); }

    private:
        Vector buff_;
        std::vector<std::shared_ptr<ForcingFunction<Vector>>> functions_;
    };

    template <class Matrix, class Vector>
    class ForcedMaterial : public ElasticMaterial<Matrix, Vector> {
    public:
        ForcedMaterial(const std::shared_ptr<ElasticMaterial<Matrix, Vector>> &material,
                       const std::shared_ptr<ForcingFunction<Vector>> &force)
            : material_(material), force_(force) {}

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            if (!material_->assemble_hessian_and_gradient(x, hessian, gradient)) {
                std::cerr << "decorated material failed to assemble" << std::endl;
                return false;
            }

            if (!force_->eval(x, force_vec_)) {
                return false;
            }

            if (material_->rescaling() != 1.0) {
                force_vec_ *= material_->rescaling();
            }

            gradient -= force_vec_;
            return true;
        }

        virtual bool stress(const Vector &x, Vector &result) override {
            if (!material_->stress(x, result)) {
                return false;
            }

            if (!force_->eval(x, force_vec_)) {
                assert(false);
                return false;
            }

            result -= force_vec_;
            return true;
        }

        bool is_linear() const override { return material_->is_linear(); }

        void clear() override { material_->clear(); }

    private:
        std::shared_ptr<ElasticMaterial<Matrix, Vector>> material_;
        std::shared_ptr<ForcingFunction<Vector>> force_;
        Vector force_vec_;
    };
}  // namespace utopia

#endif  // UTOPIA_FORCED_MATERIAL_HPP
