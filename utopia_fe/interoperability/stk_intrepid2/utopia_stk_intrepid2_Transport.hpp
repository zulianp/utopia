#ifndef UTOPIA_STK_INTREPID2_TRANSPORT_HPP
#define UTOPIA_STK_INTREPID2_TRANSPORT_HPP

#include "utopia_kokkos_Transport.hpp"
#include "utopia_stk_intrepid2_Assembler.hpp"

#include <memory>

namespace utopia {

    namespace stk {

        class Transport final : public StkIntrepid2ProxyAssembler {
        public:
            using Super = utopia::stk::StkIntrepid2ProxyAssembler;
            using Field = utopia::Field<utopia::stk::FunctionSpace>;
            using Scalar = Traits<utopia::stk::FunctionSpace>::Scalar;
            using Vector = Traits<utopia::stk::FunctionSpace>::Vector;
            using Matrix = Traits<utopia::stk::FunctionSpace>::Matrix;

            using Environment = utopia::Environment<utopia::stk::FunctionSpace>;

            Transport(const std::shared_ptr<FE> &fe);
            ~Transport();

            inline int n_vars() const override { return 1; }

            inline std::string name() const override { return "Transport"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            void read(Input &in) override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_TRANSPORT_HPP
