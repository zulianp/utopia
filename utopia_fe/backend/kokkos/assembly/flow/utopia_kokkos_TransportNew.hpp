#ifndef UTOPIA_KOKKOS_TRANSPORT_NEW_HPP
#define UTOPIA_KOKKOS_TRANSPORT_NEW_HPP

#include "utopia_fe_base.hpp"

#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_Material.hpp"

#include <memory>

namespace utopia {

    namespace kokkos {

        template <class FunctionSpace, class FE_>
        class TransportNew final : public utopia::Material<FunctionSpace, FE_> {
        public:
            using FE = FE_;
            using Super = utopia::Material<FunctionSpace, FE_>;
            using Scalar = typename Traits<FunctionSpace>::Scalar;

            inline int n_vars() const override { return 1; }
            inline std::string name() const override { return "TransportNew"; }

            inline bool has_hessian() const override { return true; }
            inline bool is_linear() const override { return true; }
            inline bool is_operator() const override { return true; };

            inline bool has_gradient() const override { return false; }
            inline bool has_value() const override { return false; }

            bool value_assemble(AssemblyMode mode) override;
            // bool gradient_assemble(AssemblyMode mode) override;
            bool hessian_assemble(AssemblyMode mode) override;

            // Matrix free hessian application
            bool apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) override;

            void read(Input &in) override;

            TransportNew();
            ~TransportNew();

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_TRANSPORT_NEW_HPP
