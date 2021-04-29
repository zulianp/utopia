#ifndef UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"
#include "utopia_intrepid2_ForwardDeclarations.hpp"

namespace utopia {
    namespace intrepid2 {

        template <class FunctionSpace>
        class OmniAssembler : public Configurable {
        public:
            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Scalar = typename Traits<FunctionSpace>::Scalar;

            OmniAssembler(const std::shared_ptr<FunctionSpace> &space);
            virtual ~OmniAssembler();

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun);
            bool assemble(const Vector &x, Matrix &jacobian);
            bool assemble(const Vector &x, Vector &fun);

            void read(Input &in) override;

            void set_environment(const std::shared_ptr<Environment<FunctionSpace>> &env);

            AssemblyMode mode() const;
            void set_mode(AssemblyMode mode);

            bool is_linear() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP
