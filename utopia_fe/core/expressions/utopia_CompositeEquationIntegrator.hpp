#ifndef UTOPIA_COMPOSITE_EQUATION_INTEGRATOR_HPP
#define UTOPIA_COMPOSITE_EQUATION_INTEGRATOR_HPP

#include "utopia_Integrator.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ModularEquationIntegrator : public EquationIntegrator<FunctionSpace> {
    public:
        using Super = utopia::EquationIntegrator<FunctionSpace>;
        using ElementMatrix = typename Super::ElementMatrix;
        using ElementVector = typename Super::ElementVector;
        using AssemblyValues = typename Super::AssemblyValues;

        using Type = ElementMatrix;
        using Super::assemble;

        static const int BACKEND_TAG = Super::BACKEND_TAG;

        using LinearIntegratorT = utopia::LinearIntegrator<FunctionSpace>;
        using BilinearIntegratorT = utopia::BilinearIntegrator<FunctionSpace>;

        virtual ~ModularEquationIntegrator() {}

        bool assemble(const Range &trial_r,
                      const Range &test_r,
                      AssemblyContext<BACKEND_TAG> &ctx,
                      ElementMatrix &mat,
                      ElementVector &vec) override {
            bool ok = bilinear_integrator_->assemble(trial_r, test_r, ctx, mat);
            assert(ok);

            if (linear_integrator_) {
                ok = linear_integrator_->assemble(test_r, ctx, vec);
                assert(ok);
            }

            return ok;
        }

        ModularEquationIntegrator(const std::shared_ptr<BilinearIntegratorT> &bilinear_integrator,
                                  const std::shared_ptr<LinearIntegratorT> &linear_integrator)
            : bilinear_integrator_(bilinear_integrator), linear_integrator_(linear_integrator) {
            assert(bilinear_integrator_);
            assert(bilinear_integrator_ || linear_integrator_);
        }

        ModularEquationIntegrator() {}

        std::shared_ptr<BilinearIntegratorT> bilinear_integrator() const { return bilinear_integrator_; }

        std::shared_ptr<LinearIntegratorT> linear_integrator() const { return linear_integrator_; }

        void bilinear_integrator(const std::shared_ptr<BilinearIntegratorT> &b) { bilinear_integrator_ = b; }

        void linear_integrator(const std::shared_ptr<LinearIntegratorT> &l) { linear_integrator_ = l; }

        int order(const AssemblyValues &ctx) const override {
            int order = 0;
            if (linear_integrator_) {
                order = std::max(linear_integrator_->order(ctx), order);
            }

            if (bilinear_integrator_) {
                order = std::max(bilinear_integrator_->order(ctx), order);
            }

            return order;
        }

        void init_values(AssemblyValues &values) const override {
            if (linear_integrator_) linear_integrator_->init_values(values);
            if (bilinear_integrator_) bilinear_integrator_->init_values(values);
        }

        const FunctionSpace &test() const override { return bilinear_integrator_->test(); }

        std::shared_ptr<FunctionSpace> test_ptr() const override { return bilinear_integrator_->test_ptr(); }

        const FunctionSpace &trial() const override { return bilinear_integrator_->trial(); }

        std::shared_ptr<FunctionSpace> trial_ptr() const override { return bilinear_integrator_->trial_ptr(); }

        inline bool is_surface() const override { return bilinear_integrator_->is_surface(); }

    private:
        std::shared_ptr<BilinearIntegratorT> bilinear_integrator_;
        std::shared_ptr<LinearIntegratorT> linear_integrator_;
    };

    template <class FS>
    std::unique_ptr<ModularEquationIntegrator<FS>> equation(const std::shared_ptr<BilinearIntegrator<FS>> &bilinear,
                                                            const std::shared_ptr<LinearIntegrator<FS>> &linear) {
        return utopia::make_unique<ModularEquationIntegrator<FS>>(bilinear, linear);
    }

    ///////////////////////////////////////////////////////////////////////////////////

    template <class FunctionSpace>
    class CompositeEquationIntegrator final : public EquationIntegrator<FunctionSpace> {
    public:
        using Super = utopia::EquationIntegrator<FunctionSpace>;
        using ElementMatrix = typename Super::ElementMatrix;
        using ElementVector = typename Super::ElementVector;
        using AssemblyValues = typename Super::AssemblyValues;

        using Type = ElementMatrix;
        using Super::assemble;

        static const int BACKEND_TAG = Super::BACKEND_TAG;

        using EquationIntegratorT = utopia::EquationIntegrator<FunctionSpace>;

        bool assemble(const Range &trial_r,
                      const Range &test_r,
                      AssemblyContext<BACKEND_TAG> &ctx,
                      ElementMatrix &mat,
                      ElementVector &vec) override {
            assert(!integrators_.empty());

            bool ok = true;
            for (auto integrator : integrators_) {
                ok = integrator->assemble(trial_r, test_r, ctx, mat, vec);
                assert(ok);
            }

            return ok;
        }

        inline void add_integrator(const std::shared_ptr<EquationIntegratorT> &integrator) {
            integrators_.push_back(integrator);
        }

        int order(const AssemblyValues &ctx) const override {
            assert(!integrators_.empty());

            int order = 0;
            for (auto integrator : integrators_) {
                order = std::max(integrator->order(ctx), order);
            }

            return order;
        }

        void init_values(AssemblyValues &values) const override {
            assert(!integrators_.empty());

            for (auto integrator : integrators_) {
                integrator->init_values(values);
            }
        }

        const FunctionSpace &test() const override {
            assert(!integrators_.empty());
            return integrators_[0]->test();
        }

        std::shared_ptr<FunctionSpace> test_ptr() const override {
            assert(!integrators_.empty());
            return integrators_[0]->test_ptr();
        }

        const FunctionSpace &trial() const override {
            assert(!integrators_.empty());
            return integrators_[0]->trial();
        }

        std::shared_ptr<FunctionSpace> trial_ptr() const override {
            assert(!integrators_.empty());
            return integrators_[0]->trial_ptr();
        }

        inline bool is_surface() const override {
            assert(!integrators_.empty());
            return integrators_[0]->is_surface();
        }

    private:
        std::vector<std::shared_ptr<EquationIntegratorT>> integrators_;
    };
}  // namespace utopia

#endif  // UTOPIA_COMPOSITE_EQUATION_INTEGRATOR_HPP
