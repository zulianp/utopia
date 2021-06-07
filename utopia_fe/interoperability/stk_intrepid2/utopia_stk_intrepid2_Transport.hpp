#ifndef UTOPIA_STK_INTREPID2_TRANSPORT_HPP
#define UTOPIA_STK_INTREPID2_TRANSPORT_HPP

#include "utopia_Field.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

#include "utopia_intrepid2_Transport.hpp"

#include "utopia_stk_FEAssembler.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include <memory>

namespace utopia {

    namespace stk {

        class StkIntrepid2Assembler
            : public utopia::intrepid2::FEAssembler<typename Traits<stk::FunctionSpace>::Scalar> {
        public:
            using Super = utopia::intrepid2::FEAssembler<typename Traits<stk::FunctionSpace>::Scalar>;
            //     using Super = utopia::FEAssembler<utopia::stk::FunctionSpace>;
            //     using Field = utopia::Field<utopia::stk::FunctionSpace>;
            //     using Scalar = Traits<utopia::stk::FunctionSpace>::Scalar;
            //     using Vector = Traits<utopia::stk::FunctionSpace>::Vector;
            //     using Matrix = Traits<utopia::stk::FunctionSpace>::Matrix;

            using Environment = utopia::Environment<utopia::stk::FunctionSpace>;

            StkIntrepid2Assembler(const std::shared_ptr<FE> &fe);
            virtual ~StkIntrepid2Assembler();

            void set_environment(const std::shared_ptr<Environment> &env);
            std::shared_ptr<Environment> environment() const;

            void set_space(const std::shared_ptr<FunctionSpace> &space);
            std::shared_ptr<FunctionSpace> space() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

        class Transport final : public StkIntrepid2Assembler {
        public:
            using Super = utopia::stk::StkIntrepid2Assembler;
            using Field = utopia::Field<utopia::stk::FunctionSpace>;
            using Scalar = Traits<utopia::stk::FunctionSpace>::Scalar;
            using Vector = Traits<utopia::stk::FunctionSpace>::Vector;
            using Matrix = Traits<utopia::stk::FunctionSpace>::Matrix;

            using Environment = utopia::Environment<utopia::stk::FunctionSpace>;

            using Intrepid2Assembler = utopia::intrepid2::FEAssembler<Scalar>;
            using Intrepid2FE = intrepid2::FE<Scalar>;
            using TensorAccumulator = Intrepid2Assembler::TensorAccumulator;

            Transport(const std::shared_ptr<FE> &fe);
            ~Transport();

            inline int n_vars() const override { return 1; }

            inline std::string name() const override { return "Transport"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            void read(Input &in) override;

            bool apply(const DynRankView &x, DynRankView &y) override;
            bool assemble_matrix() override;

            void set_matrix_accumulator(const std::shared_ptr<TensorAccumulator> &matrix_accumulator) override;

            void set_vector_accumulator(const std::shared_ptr<TensorAccumulator> &vector_accumulator) override;

            void set_scalar_accumulator(const std::shared_ptr<TensorAccumulator> &scalar_accumulator) override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_TRANSPORT_HPP
