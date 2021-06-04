#ifndef UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP

// Utopia
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

// UtopiaFE
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

// Utopia::Stk
#include "utopia_stk_FEAssembler.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"
#include "utopia_stk_FunctionSpace.hpp"

namespace utopia {

    template <>
    class OmniAssembler<utopia::stk::FunctionSpace> final : public FEAssembler<utopia::stk::FunctionSpace> {
    public:
        using FunctionSpace = utopia::stk::FunctionSpace;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Environment = utopia::Environment<FunctionSpace>;

        ~OmniAssembler();
        OmniAssembler(const std::shared_ptr<FunctionSpace> &space);

        void clear() override;
        std::string name() const override;
        bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override;

        bool assemble(const Vector &x, Matrix &hessian) override;
        bool assemble(const Vector &x, Vector &gradient) override;

        // For linear only
        bool assemble(Matrix &hessian) override;
        bool assemble(Vector &gradient) override;

        void read(Input &in) override;

        void set_environment(const std::shared_ptr<Environment> &env) override;
        std::shared_ptr<Environment> environment() const override;

        void set_space(const std::shared_ptr<FunctionSpace> &space) override;
        std::shared_ptr<FunctionSpace> space() const override;

        bool is_linear() const override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP
