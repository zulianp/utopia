#ifndef UTOPIA_UNUSPPORTED_ALGORITHMS_HPP
#define UTOPIA_UNUSPPORTED_ALGORITHMS_HPP

#include "utopia_IsNotSupported.hpp"
#include "utopia_libmesh_ForwardDeclarations.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ImplicitObstacle;

    template <class FunctionSpace>
    class AnalyticObstacle;

    template <>
    class IsNotSupported<ImplicitObstacle<libmesh::FunctionSpace>> {
    public:
        static constexpr bool value{true};
    };

    template <>
    class IsNotSupported<AnalyticObstacle<libmesh::FunctionSpace>> {
    public:
        static constexpr bool value{true};
    };

}  // namespace utopia

#endif  // UTOPIA_UNUSPPORTED_ALGORITHMS_HPP
