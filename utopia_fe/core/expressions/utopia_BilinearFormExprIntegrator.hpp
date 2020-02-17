#ifndef UTOPIA_BILINEAR_FORM_EXPR_INTEGRATOR_HPP
#define UTOPIA_BILINEAR_FORM_EXPR_INTEGRATOR_HPP

#include "utopia_Integrator.hpp"
#include "utopia_FEEval_Integrator.hpp"

namespace utopia {

    template<class FunctionSpace, class Expr>
    class BiLinearFormExprIntegrator final : public BilinearIntegrator<FunctionSpace> {
    public:
        using Super = utopia::BilinearIntegrator<FunctionSpace>;
        static const int BACKEND_TAG = Super::BACKEND_TAG;
        using ElementMatrix  = typename Super::ElementMatrix;
        using AssemblyValues = typename Super::AssemblyValues;
        using Super::assemble;

        BiLinearFormExprIntegrator(const Expr &expr)
        : expr_(expr), space_(find_test_space<FunctionSpace>(expr))
        {
            assert(space_);
        }

        BiLinearFormExprIntegrator(const std::shared_ptr<FunctionSpace> &space, const Expr &expr)
        : expr_(expr), space_(space)
        {
            assert(space_);
        }

        bool assemble(
            const Range &trial_r,
            const Range &test_r,
            AssemblyContext<BACKEND_TAG> &ctx,
            ElementMatrix &result) override
         {
            UTOPIA_UNUSED(trial_r); //handled internally
            UTOPIA_UNUSED(test_r); //handled internally
            result += utopia::eval(expr_, ctx);
            return true;
         }

        int order(const AssemblyValues &ctx) const override
        {
            return FunctionalTraits<Expr, AssemblyValues>::order(expr_, ctx);
        }

        void init_values(AssemblyValues &values) const override
        {
            values.init_fe_flags(expr_);
        }

        const FunctionSpace &test() const override
        {
            assert(space_);
            return *space_;
        }

        std::shared_ptr<FunctionSpace> test_ptr() const override
        {
            return space_;
        }

        bool is_surface() const override
        {
            return expr_.is_surface();
        }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
        std::shared_ptr<FunctionSpace> space_;
    };

    template<class FunctionSpace, class Expr>
    std::unique_ptr<BilinearIntegrator<FunctionSpace>> bilinear_form(
        const std::shared_ptr<FunctionSpace> &fs,
        const Expr &expr
        )
    {
        return utopia::make_unique< BiLinearFormExprIntegrator<FunctionSpace, Expr> >(
            fs,
            expr
        );
    }

    template<class FS, class Expr>
    std::unique_ptr<BilinearIntegrator<FS>> bilinear_form(
        FunctionSpace<FS> &fs,
        const Expr &expr
        )
    {
        return utopia::make_unique< BiLinearFormExprIntegrator<FS, Expr> >(
            make_ref(fs.derived()),
            expr
        );
    }

}

#endif //UTOPIA_BILINEAR_FORM_EXPR_INTEGRATOR_HPP
