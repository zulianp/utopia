#ifndef UTOPIA_STENCIL_VIEW_HPP
#define UTOPIA_STENCIL_VIEW_HPP

#include "utopia_MeshView.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    template <class Mesh>
    class StencilView {};

    template <int Dim>
    class TensorShift {
    public:
    };

    template <>
    class TensorShift<2> {
    public:
        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static bool apply(const Dims &dims, const int &idx, Array &array) {
            // sorted indexing
            switch (idx) {
                case 0: {
                    if (array[0] == 0) return false;
                    array[0] -= 1;
                    return true;
                }
                case 1: {
                    if (array[1] == 0) return false;
                    array[1] -= 1;
                    return true;
                }
                case 2: {
                    // identity
                    return true;
                }
                case 3: {
                    if (array[0] == dims[0]) return false;
                    array[0] += 1;
                    return true;
                }
                case 4: {
                    if (array[1] == dims[1]) return false;
                    array[1] += 1;
                    return true;
                }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return false;
                }
            }
        }
    };

    template <class Elem, class ExecutionSpace, typename... Args>
    class StencilView<TensorMeshView<Elem, ExecutionSpace, Args...>> {
    public:
        static const int Dim = Elem::Dim;

        using Mesh = utopia::TensorMeshView<Elem, ExecutionSpace, Args...>;
        using SizeType = typename Mesh::SizeType;
        using TensorShift = utopia::TensorShift<Dim>;
        UTOPIA_INLINE_FUNCTION constexpr static SizeType invalid() { return -1; }

        StencilView(const Mesh &mesh) : mesh_(mesh) {}

        UTOPIA_INLINE_FUNCTION constexpr static SizeType size() { return Dim * 2 + 1; }

        UTOPIA_INLINE_FUNCTION SizeType index(const SizeType &i, const SizeType &k) const {
            ArrayView<SizeType, Dim> tensor_index;
            mesh_.node_linear_to_tensor_index(i, tensor_index);
            if (!TensorShift::apply(mesh_.dims(), k, tensor_index)) {
                return invalid();
            }

            return mesh_.node_tensor_to_linear_index(tensor_index);
        }

    private:
        Mesh mesh_;
    };

}  // namespace utopia

#endif  // UTOPIA_STENCIL_VIEW_HPP
