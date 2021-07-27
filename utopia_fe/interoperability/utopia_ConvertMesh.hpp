#ifndef UTOPIA_CONVERT_MESH_HPP
#define UTOPIA_CONVERT_MESH_HPP

namespace utopia {
    template <class MeshIn, class MeshOut>
    class ConvertMesh {};

    // FIXME use convert once we have static polymorphism up
    template <class MeshIn, class MeshOut>
    inline void convert_mesh(const MeshIn &in, MeshOut &out) {
        ConvertMesh<MeshIn, MeshOut>::apply(in, out);
    }
}  // namespace utopia

#endif  // UTOPIA_CONVERT_MESH_HPP
