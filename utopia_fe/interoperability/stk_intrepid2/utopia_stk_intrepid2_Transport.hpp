#ifndef UTOPIA_STK_INTREPID2_TRANSPORT_HPP
#define UTOPIA_STK_INTREPID2_TRANSPORT_HPP

#include "utopia_Field.hpp"
#include "utopia_fe_Environment.hpp"
#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_stk_FEAssembler.hpp"

#include <memory>

namespace utopia {
    namespace stk {

        class StkIntrepid2Assembler : public FEAssembler<utopia::stk::FunctionSpace> {
        public:
            using Super = utopia::FEAssembler<utopia::stk::FunctionSpace>;
            using Field = utopia::Field<utopia::stk::FunctionSpace>;
            using Scalar = Traits<utopia::stk::FunctionSpace>::Scalar;
            using Intrepid2Assembler = utopia::intrepid2::FEAssembler<Scalar>;
            using Intrepid2FE = intrepid2::FE<Scalar>;
            using TensorAccumulator = Intrepid2Assembler::TensorAccumulator;

            virtual ~StkIntrepid2Assembler();
            StkIntrepid2Assembler();

            bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) final;

            void set_accumulator(const std::shared_ptr<TensorAccumulator> &accumulator);
            inline std::shared_ptr<TensorAccumulator> accumulator();

            void set_assembler(const std::shared_ptr<Intrepid2Assembler> &assembler);
            std::shared_ptr<Intrepid2Assembler> assembler();

            void set_fe(const std::shared_ptr<Intrepid2FE> &fe);
            std::shared_ptr<Intrepid2FE> fe();
            void ensure_fe(const int quadrature_order);

            void read(Input &in) override;

            virtual bool assemble_element_tensors();
            inline AssemblyMode assembly_mode() const { return mode_; }
            inline void set_assembly_mode(AssemblyMode mode) { mode_ = mode; }

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
            AssemblyMode mode_{ADD_MODE};
        };

        class Transport : public StkIntrepid2Assembler {
        public:
            using Super = utopia::stk::StkIntrepid2Assembler;

            ~Transport();
            Transport();

            inline std::string name() const override { return "Transport"; }
            void read(Input &in) override;

            // bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

        class Mass : public StkIntrepid2Assembler {
        public:
            using Super = utopia::stk::StkIntrepid2Assembler;

            ~Mass();
            Mass();

            inline std::string name() const override { return " Mass"; }
            void read(Input &in) override;

            // bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override;

            bool assemble_element_tensors() override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            void init();
            void ensure_assembler();
        };

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_TRANSPORT_HPP
