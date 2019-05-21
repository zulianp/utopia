#ifndef UTOPIA_COMPOSITE_EQUATION_INTEGRATOR_HPP
#define UTOPIA_COMPOSITE_EQUATION_INTEGRATOR_HPP

#include "utopia_Integrator.hpp"

namespace utopia {

    template<class FunctionSpace>
    class ModularEquationIntegrator final : public EquationIntegrator<FunctionSpace> {
    public:
        using Super = utopia::EquationIntegrator<FunctionSpace>;
        using ElementMatrix  = typename Super::ElementMatrix;
        using ElementVector  = typename Super::ElementVector;
        using AssemblyValues = typename Super::AssemblyValues;

        using Type = ElementMatrix;
        using Super::assemble;

        static const int BACKEND_TAG = Super::BACKEND_TAG;

        using LinearIntegratorT   = utopia::LinearIntegrator<FunctionSpace>;
        using BilinearIntegratorT = utopia::BilinearIntegrator<FunctionSpace>;

        virtual ~ModularEquationIntegrator() {}
        
        virtual bool assemble(
            const Range &trial_r,
            const Range &test_r,
            AssemblyContext<BACKEND_TAG> &ctx,
            ElementMatrix &mat,
            ElementVector &vec) override
        {
            bool ok = bilinear_integrator_->assemble(trial_r, test_r, ctx, mat); assert(ok);
                 ok = linear_integrator_->assemble(test_r, ctx, vec);            assert(ok);

            return ok;
        }

        ModularEquationIntegrator(
            const std::shared_ptr<BilinearIntegratorT> &bilinear_integrator,
            const std::shared_ptr<LinearIntegratorT>   &linear_integrator)
        : bilinear_integrator_(bilinear_integrator),
          linear_integrator_(linear_integrator)
        {
            assert(bilinear_integrator_);
            assert(linear_integrator_);
        }

        int order(const AssemblyValues &ctx) const override
        {
            return std::max(linear_integrator_->order(ctx), bilinear_integrator_->order(ctx));
        }

        void init_values(AssemblyValues &values) const override
        {
            linear_integrator_->init_values(values);
            bilinear_integrator_->init_values(values);
        }

        const FunctionSpace &test() const override
        {
            return bilinear_integrator_->test();
        }

        std::shared_ptr<FunctionSpace> test_ptr() const override
        {
            return bilinear_integrator_->test_ptr();
        }

        const FunctionSpace &trial() const override
        {
            return bilinear_integrator_->trial();
        }

        std::shared_ptr<FunctionSpace> trial_ptr() const override
        {
            return bilinear_integrator_->trial_ptr();
        }

        inline bool is_surface() const override
        {
            return bilinear_integrator_->is_surface() || linear_integrator_->is_surface(); 
        }

    private:
        std::shared_ptr<LinearIntegratorT>   linear_integrator_;
        std::shared_ptr<BilinearIntegratorT> bilinear_integrator_;
    };


    template<class FS>
    std::unique_ptr<ModularEquationIntegrator<FS>> equation(
        const std::shared_ptr<BilinearIntegrator<FS>> &bilinear,
        const std::shared_ptr<LinearIntegrator<FS>> &linear)
    {
        return utopia::make_unique<ModularEquationIntegrator<FS>>(bilinear, linear);
    }

}

#endif //UTOPIA_COMPOSITE_EQUATION_INTEGRATOR_HPP
