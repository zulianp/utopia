#ifndef UTOPIA_KOKKOS_OMNIMATERIAL_HPP
#define UTOPIA_KOKKOS_OMNIMATERIAL_HPP

#include "utopia_Input.hpp"
#include "utopia_Material.hpp"
#include "utopia_SimulationTime.hpp"
#include "utopia_Traits.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"
#include "utopia_kokkos_Field.hpp"

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE>
        class OmniMaterial final : public utopia::Material<FunctionSpace, FE> {
        public:
            using Super = utopia::Material<FunctionSpace, FE>;
            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Scalar = typename Traits<FunctionSpace>::Scalar;
            using Environment = typename Traits<FunctionSpace>::Environment;
            using SimulationTime = utopia::SimulationTime<Scalar>;

            OmniMaterial(const std::shared_ptr<FunctionSpace> &space);
            ~OmniMaterial();

            void read(Input &in) override;

            int n_vars() const override;
            std::string name() const override;

            bool has_hessian() const override;
            bool is_linear() const override;
            bool is_operator() const override;

            bool has_gradient() const override;
            bool has_value() const override;

            bool value_assemble(AssemblyMode mode) override;
            bool gradient_assemble(AssemblyMode mode) override;
            bool hessian_assemble(AssemblyMode mode) override;
            bool apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) override;

            void set_time(const std::shared_ptr<SimulationTime> &time) override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_OMNIMATERIAL_HPP
