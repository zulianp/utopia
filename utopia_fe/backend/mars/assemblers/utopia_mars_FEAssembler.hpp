#ifndef UTOPIA_MARS_FE_ASSEMBLER_HPP
#define UTOPIA_MARS_FE_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_FEAssembler.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_mars_ForwardDeclarations.hpp"
#include "utopia_mars_FunctionSpace.hpp"

namespace utopia {
    namespace mars {

        class FEAssembler : public Describable, public utopia::FEAssembler<utopia::mars::FunctionSpace> {
        public:
            using Super = utopia::FEAssembler<utopia::mars::FunctionSpace>;
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;

            virtual ~FEAssembler() = default;
            virtual bool assemble(const Vector &x, Matrix &mat, Vector &vec) = 0;
            virtual bool assemble(const Vector &x, Matrix &mat) = 0;
            virtual bool assemble(const Vector &x, Vector &vec) = 0;
            virtual bool assemble(Matrix &mat) = 0;
            virtual void set_environment(const std::shared_ptr<Environment> &env) = 0;
            virtual bool is_linear() const = 0;

            // TODO
            bool is_operator() const override { return false; }

            inline void set_mode(AssemblyMode mode) { mode_ = mode; }
            inline AssemblyMode mode() const { return mode_; }

            virtual void set_space(const std::shared_ptr<FunctionSpace> &space) = 0;

        private:
            AssemblyMode mode_{ADD_MODE};
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_FE_ASSEMBLER_HPP
