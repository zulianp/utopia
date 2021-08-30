#ifndef UTOPIA_STK_INTREPID2_ASSEMBLER_HPP
#define UTOPIA_STK_INTREPID2_ASSEMBLER_HPP

#include "utopia_Field.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Transport.hpp"

#include "utopia_intrepid2_FE.hpp"

#include "utopia_stk_FEAssembler.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include <memory>

namespace utopia {

    namespace stk {

        class StkIntrepid2Assembler
            : public utopia::kokkos::FEAssembler<intrepid2::FE<Traits<stk::FunctionSpace>::Scalar>> {
        public:
            using Scalar = Traits<stk::FunctionSpace>::Scalar;
            using Super = utopia::kokkos::FEAssembler<intrepid2::FE<Traits<stk::FunctionSpace>::Scalar>>;
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

        class StkIntrepid2ProxyAssembler : public StkIntrepid2Assembler {
        public:
            using Super = utopia::stk::StkIntrepid2Assembler;
            using Intrepid2Assembler = utopia::kokkos::FEAssembler<intrepid2::FE<Traits<stk::FunctionSpace>::Scalar>>;
            using Intrepid2FE = intrepid2::FE<Scalar>;
            using MatrixAccumulator = Intrepid2Assembler::MatrixAccumulator;
            using VectorAccumulator = Intrepid2Assembler::VectorAccumulator;
            using ScalarAccumulator = Intrepid2Assembler::ScalarAccumulator;

            using VectorView = Intrepid2Assembler::VectorView;

            virtual ~StkIntrepid2ProxyAssembler();
            StkIntrepid2ProxyAssembler(const std::shared_ptr<FE> &fe);

            void set_matrix_accumulator(const std::shared_ptr<MatrixAccumulator> &matrix_accumulator) override;

            void set_vector_accumulator(const std::shared_ptr<VectorAccumulator> &vector_accumulator) override;

            void set_scalar_accumulator(const std::shared_ptr<ScalarAccumulator> &scalar_accumulator) override;

            void ensure_matrix_accumulator() override;
            void ensure_vector_accumulator() override;
            void ensure_scalar_accumulator() override;

            bool apply(const VectorView &x, VectorView &y) override;
            bool assemble_matrix() override;

            void set_assembler(const std::shared_ptr<Intrepid2Assembler> &assembler);
            const std::shared_ptr<Intrepid2Assembler> &assembler() const;

            void set_time(const std::shared_ptr<SimulationTime> &time);

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_ASSEMBLER_HPP
