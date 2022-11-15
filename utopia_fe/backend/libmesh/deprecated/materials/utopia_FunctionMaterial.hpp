#ifndef UTOPIA_FUNCTION_MATERIAL_HPP
#define UTOPIA_FUNCTION_MATERIAL_HPP

#include "utopia_ElasticMaterial.hpp"
#include "utopia_Function.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class FunctionMaterial : public ElasticMaterial<Matrix, Vector> {
    public:
        FunctionMaterial(const std::shared_ptr<Function<Matrix, Vector> > &fun) : fun_(fun) {}

        virtual ~FunctionMaterial() {}

        virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) {
            bool ok = fun_->gradient(x, gradient);
            ok &= fun_->hessian(x, hessian);
            return ok;
        }

        bool assemble_hessian(const Vector &x, Matrix &H) override { return fun_->hessian(x, H); }
        bool assemble_gradient(const Vector &, Vector &g) override { return fun_->gradient(x, g); }

        bool assemble_hessian(Matrix &) override {
            assert(false);
            return false;
        }

        virtual void clear() {}

    private:
        std::shared_ptr<Function<Matrix, Vector> > fun_;
    };
}  // namespace utopia

#endif  // UTOPIA_FUNCTION_MATERIAL_HPP
