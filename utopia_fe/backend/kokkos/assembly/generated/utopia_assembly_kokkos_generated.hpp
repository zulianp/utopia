#ifndef UTOPIA_ASSEMBLY_KOKKOS_GENERATED_HPP
#define UTOPIA_ASSEMBLY_KOKKOS_GENERATED_HPP

namespace utopia {
    namespace kernels {

        // Forward declaratins
        template <class Elem>
        class Mass;

        template <class Elem>
        class LplaceOperator;

    }  // namespace kernels

    namespace kokkos {

        template <class FunctionSpace, class FE>
        class AssemblerRegistry;

        template <class FunctionSpace, class FE>
        void register_generated_assemblers(AssemblerRegistry<FunctionSpace, FE> &registry);

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_ASSEMBLY_KOKKOS_GENERATED_HPP
