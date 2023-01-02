#ifndef UTOPIA_KOKKOS_OMNIMATERIAL_IMPL_HPP
#define UTOPIA_KOKKOS_OMNIMATERIAL_IMPL_HPP

#include "utopia_kokkos_OmniMaterial.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE>
        class OmniMaterial<FunctionSpace, FE>::Impl {};

        template <class FunctionSpace, class FE>
        OmniMaterial<FunctionSpace, FE>::OmniMaterial(const std::shared_ptr<FunctionSpace> &space)
            : impl_(utopia::make_unique<Impl>()) {}

        template <class FunctionSpace, class FE>
        OmniMaterial<FunctionSpace, FE>::~OmniMaterial() {}

        template <class FunctionSpace, class FE>
        void OmniMaterial<FunctionSpace, FE>::read(Input &in) {}

        template <class FunctionSpace, class FE>
        int OmniMaterial<FunctionSpace, FE>::n_vars() const {}

        template <class FunctionSpace, class FE>
        std::string OmniMaterial<FunctionSpace, FE>::name() const {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::has_hessian() const {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::is_linear() const {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::is_operator() const {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::has_gradient() const {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::has_value() const {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::value_assemble(AssemblyMode mode) {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::gradient_assemble(AssemblyMode mode) {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::hessian_assemble(AssemblyMode mode) {}

        template <class FunctionSpace, class FE>
        bool OmniMaterial<FunctionSpace, FE>::apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) {}

        template <class FunctionSpace, class FE>
        void OmniMaterial<FunctionSpace, FE>::set_time(const std::shared_ptr<SimulationTime> &time) {}

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_OMNIMATERIAL_IMPL_HPP
