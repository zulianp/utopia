#ifndef UTOPIA_KOKKOS_MATERIAL_FACTORY_HPP
#define UTOPIA_KOKKOS_MATERIAL_FACTORY_HPP

#include "utopia_Input.hpp"
#include "utopia_Material.hpp"

#include <memory>

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE>
        class MaterialFactory {
        public:
            static std::unique_ptr<utopia::Material<FunctionSpace, FE>> make(const int ndims, const std::string &name);
            class Impl;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MATERIAL_FACTORY_HPP
