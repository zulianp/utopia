#ifndef UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_FEAssembler.hpp"

#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"
#include "utopia_intrepid2_ForwardDeclarations.hpp"

namespace utopia {
    namespace intrepid2 {

        template <class FunctionSpace>
        class OmniAssembler : public FEAssembler<FunctionSpace> {
        public:
            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Scalar = typename Traits<FunctionSpace>::Scalar;

            using Environment = utopia::Environment<FunctionSpace>;

            OmniAssembler(const std::shared_ptr<FunctionSpace> &space);
            virtual ~OmniAssembler();

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun) override;
            bool assemble(const Vector &x, Matrix &jacobian) override;
            bool assemble(const Vector &x, Vector &fun) override;

            // For linear only
            bool assemble(Matrix &jacobian) override;
            bool assemble(Vector &fun) override;

            void read(Input &in) override;

            AssemblyMode mode() const;
            void set_mode(AssemblyMode mode);

            std::string name() const override;

            void set_environment(const std::shared_ptr<Environment> &env) override;
            std::shared_ptr<Environment> environment() const override;
            void set_space(const std::shared_ptr<FunctionSpace> &space) override;
            std::shared_ptr<FunctionSpace> space() const override;
            bool is_linear() const override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_OMNI_ASSEMBLER_HPP
