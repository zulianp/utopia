#ifndef UTOPIA_KOKKOS_FE_FORWARD_DECLARATIONS_HPP
#define UTOPIA_KOKKOS_FE_FORWARD_DECLARATIONS_HPP

namespace utopia {
    namespace kokkos {
        template <typename Scalar_, class Traits>
        class FE;

        template <typename FE_>
        class Field;
    }  // namespace kokkos

    template <typename T, class FETraits>
    class Traits<utopia::kokkos::FE<T, FETraits>> {
    public:
        using Scalar = T;
        using Field = utopia::kokkos::Field<utopia::kokkos::FE<T, FETraits>>;
    };
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_FE_FORWARD_DECLARATIONS_HPP
