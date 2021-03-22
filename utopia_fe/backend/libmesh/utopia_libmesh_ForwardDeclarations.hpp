#ifndef UTOPIA_LIBMESH_FORWARD_DECLARATIONS_HPP
#define UTOPIA_LIBMESH_FORWARD_DECLARATIONS_HPP

namespace libMesh {
    class MeshBase;
    class EquationSystems;
    class DofMap;

    template <typename T>
    class NumericVector;

    namespace Parallel {
        class Communicator;
    }

}  // namespace libMesh

namespace utopia {
    namespace libmesh {
        class Mesh;
        class MeshInitializer;
        class FunctionSpace;

        class FunctionSubspace;
    }  // namespace libmesh
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_FORWARD_DECLARATIONS_HPP
