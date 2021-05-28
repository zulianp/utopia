#ifndef UTOPIA_MARS_OMNIASSEMBLER_HPP
#define UTOPIA_MARS_OMNIASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_mars_ForwardDeclarations.hpp"
#include "utopia_mars_FunctionSpace.hpp"

namespace utopia {
    namespace mars {

        class OmniAssembler : public Configurable {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;

            OmniAssembler(const std::shared_ptr<mars::FunctionSpace> &space);
            virtual ~OmniAssembler();

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun);
            bool assemble(const Vector &x, Matrix &jacobian);
            bool assemble(const Vector &x, Vector &fun);

            bool assemble(Matrix &jacobian);
            void read(Input &in) override;

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &env);

            bool is_linear() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace mars

    template <>
    class OmniAssembler<utopia::mars::FunctionSpace> final : public utopia::mars::OmniAssembler {
    public:
        using Super = utopia::mars::OmniAssembler;
        using Super::Super;
    };

}  // namespace utopia

#endif  // UTOPIA_MARS_OMNIASSEMBLER_HPP