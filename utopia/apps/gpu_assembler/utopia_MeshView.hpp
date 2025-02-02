#ifndef UTOPIA_MESH_VIEW_HPP
#define UTOPIA_MESH_VIEW_HPP

#include "utopia_Base.hpp"
#include "utopia_Box.hpp"
#include "utopia_MemType.hpp"
#include "utopia_StaticMath.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Views.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <array>

namespace utopia {

    template <int Dim>
    class CellToNodeIndex {};

    template <class Scalar, typename... Args>
    class Box<Kokkos::DualView<Scalar *, Args...>> {
    public:
        using View = Kokkos::DualView<Scalar *, Args...>;

        Box(const int dim) : min("box_min", dim), max("box_max", dim) {}

        template <std::size_t UDim>
        explicit Box(const Box<std::array<Scalar, UDim>> &box) : min("box_min", UDim), max("box_max", UDim) {
            auto h_min = min.view_host();
            auto h_max = max.view_host();

            static const int Dim = UDim;

            for (int i = 0; i < Dim; ++i) {
                h_min[i] = box.min[i];
                h_max[i] = box.max[i];
            }
        }

        View min;
        View max;
    };

    template <>
    class CellToNodeIndex<2> {
    public:
        using SizeType = std::size_t;

        template <class Dims, class TensorIndexView, class NodeIndexView>
        UTOPIA_INLINE_FUNCTION static void build(const Dims &dims,
                                                 const TensorIndexView &cell_tensor_idx,
                                                 NodeIndexView &node_idx) {
            const SizeType i = cell_tensor_idx[0];
            const SizeType j = cell_tensor_idx[1];
            const SizeType n_p = dims[1] + 1;

            // flipped for dof-consistency and ccw local

            node_idx[0] = i * n_p + j;
            node_idx[3] = i * n_p + (j + 1);
            node_idx[2] = (i + 1) * n_p + (j + 1);
            node_idx[1] = (i + 1) * n_p + j;
        }
    };

    template <class Mesh>
    class TensorMeshImpl {
    public:
        static const int Dim = Mesh::Dim;

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType local_element_index_begin(const Dims &dims,
                                                                         const Array &local_elements_begin) {
            return element_tensor_to_linear_index(dims, local_elements_begin);
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType local_element_index_end(const Dims &dims,
                                                                       const Array &local_elements_end) {
            SizeType result = local_elements_end[Dim - 1] - 1;

            for (int i = Dim - 2; i >= 0; --i) {
                result *= dims[i];
                result += local_elements_end[i] - 1;
            }

            return result + 1;
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static void element_linear_to_tensor_index(const Dims &dims,
                                                                          const SizeType &element_idx,
                                                                          Array &tensor_index) {
            SizeType current = element_idx;

            for (int i = 0; i < Dim; ++i) {
                const int next = current / dims[i];
                tensor_index[i] = current - next * dims[i];
                current = next;
            }

            // UTOPIA_DEVICE_ASSERT(element_tensor_to_linear_index(tensor_index) == element_idx);
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType element_tensor_to_linear_index(const Dims &dims,
                                                                              const Array &tensor_index) {
            SizeType result = tensor_index[Dim - 1];

            for (int i = Dim - 2; i >= 0; --i) {
                result *= dims[i];
                result += tensor_index[i];
            }

            return result;
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType local_node_index_begin(const Dims &dims,
                                                                      const Array &local_elements_begin) {
            return node_tensor_to_linear_index(dims, local_elements_begin);
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType local_node_index_end(const Dims &dims, const Array &local_elements_end) {
            SizeType result = local_elements_end[Dim - 1] - 1 + (local_elements_end[Dim - 1] == dims[Dim - 1]);

            for (int i = Dim - 2; i >= 0; --i) {
                result *= (dims[i] + 1);
                result += local_elements_end[i] - 1 + (local_elements_end[i] == dims[i]);
            }

            return result + 1;
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static void node_linear_to_tensor_index(const Dims &dims,
                                                                       const SizeType &node_idx,
                                                                       Array &tensor_index) {
            SizeType current = node_idx;

            for (int i = 0; i < Dim; ++i) {
                const SizeType next = current / (dims[i] + 1);
                tensor_index[i] = current - next * (dims[i] + 1);
                current = next;
            }

            // UTOPIA_DEVICE_ASSERT(node_tensor_to_linear_index(tensor_index) == node_idx);
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType node_tensor_to_linear_index(const Dims &dims,
                                                                           const Array &tensor_index) {
            SizeType result = tensor_index[0];

            for (SizeType i = 1; i < Dim; ++i) {
                result *= dims[i] + 1;
                result += tensor_index[i];
            }

            return result;
        }
    };

    template <class Mesh>
    class TensorMeshImplXMajor {
    public:
        static const int Dim = Mesh::Dim;

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType local_element_index_begin(const Dims &dims,
                                                                         const Array &local_elements_begin) {
            return element_tensor_to_linear_index(dims, local_elements_begin);
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType local_element_index_end(const Dims &dims,
                                                                       const Array &local_elements_end) {
            SizeType result = local_elements_end[0] - 1;

            for (int i = 1; i < Dim; ++i) {
                result *= dims[i];
                result += local_elements_end[i] - 1;
            }

            return result + 1;
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static void element_linear_to_tensor_index(const Dims &dims,
                                                                          const SizeType &element_idx,
                                                                          Array &tensor_index) {
            SizeType current = element_idx;
            const int last = Dim - 1;

            for (int i = last; i >= 0; --i) {
                const int next = current / dims[i];
                tensor_index[i] = current - next * dims[i];
                current = next;
            }

            UTOPIA_DEVICE_ASSERT(element_tensor_to_linear_index(tensor_index) == element_idx);
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType element_tensor_to_linear_index(const Dims &dims,
                                                                              const Array &tensor_index) {
            SizeType result = tensor_index[0];

            for (int i = 1; i < Dim; ++i) {
                result *= dims[i];
                result += tensor_index[i];
            }

            return result;
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType local_node_index_begin(const Dims &dims,
                                                                      const Array &local_elements_begin) {
            return node_tensor_to_linear_index(dims, local_elements_begin);
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType local_node_index_end(const Dims &dims, const Array &local_elements_end) {
            SizeType result = local_elements_end[0] - 1 + (local_elements_end[0] == dims[0]);

            for (int i = 1; i < Dim; ++i) {
                result *= (dims[i] + 1);
                result += local_elements_end[i] - 1 + (local_elements_end[i] == dims[i]);
            }

            return result + 1;
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static void node_linear_to_tensor_index(const Dims &dims,
                                                                       const SizeType &node_idx,
                                                                       Array &tensor_index) {
            SizeType current = node_idx;
            const int last = Dim - 1;

            for (int i = last; i >= 0; --i) {
                const SizeType next = current / (dims[i] + 1);
                tensor_index[i] = current - next * (dims[i] + 1);
                current = next;
            }

            UTOPIA_DEVICE_ASSERT(node_tensor_to_linear_index(tensor_index) == node_idx);
        }

        template <class Dims, class Array>
        UTOPIA_INLINE_FUNCTION static SizeType node_tensor_to_linear_index(const Dims &dims,
                                                                           const Array &tensor_index) {
            SizeType result = tensor_index[0];

            for (SizeType i = 1; i < Dim; ++i) {
                result *= dims[i] + 1;
                result += tensor_index[i];
            }

            return result;
        }
    };

    template <class Elem_, class ExecutionSpace, typename... Args>
    class TensorMeshView {
    public:
        using SizeType = std::size_t;
        using Elem = Elem_;

        static const int Dim = Elem::Dim;
        static const int NNodesXElem = utopia::Power<Dim, 2>::value;

        static const int LEFT = -1;
        static const int RIGHT = 1;

        using Scalar = typename Elem::Scalar;
        using NodeIndexView = utopia::ArrayView<SizeType, NNodesXElem>;
        using TensorIndexView = utopia::ArrayView<SizeType, Dim>;
        using IndexView = Kokkos::View<SizeType *, ExecutionSpace>;
        using ScalarView = Kokkos::View<Scalar *, ExecutionSpace>;

        template <class Array>
        UTOPIA_INLINE_FUNCTION void nodes(const SizeType &element_idx, Array &node_idx) const {
            TensorIndexView tensor_index;
            element_linear_to_tensor_index(element_idx, tensor_index);
            CellToNodeIndex<Dim>::build(dims_, tensor_index, node_idx);
        }

        template <class Point>
        UTOPIA_INLINE_FUNCTION void point(const SizeType &point_idx, Point &point) {
            TensorIndexView tensor_index;
            node_linear_to_tensor_index(point_idx, tensor_index);

            for (int d = 0; d < Dim; ++d) {
                point[d] = tensor_index[d] * h_[d] + box_.min[d];
            }
        }

        UTOPIA_INLINE_FUNCTION void elem(const SizeType &element_idx, Elem &elem) const {
            TensorIndexView tensor_index;
            element_linear_to_tensor_index(element_idx, tensor_index);

            StaticVector2<Scalar> translation;
            for (int d = 0; d < Dim; ++d) {
                translation(d) = tensor_index[d] * h_[d] + box_.min[d];
            }

            elem.set(translation, h_);
            nodes(element_idx, elem.nodes());
        }

        TensorMeshView(const IndexView &dims,
                       const IndexView &local_elements_begin,
                       const IndexView &local_elements_end,
                       const Box<ScalarView> &box,
                       const ScalarView &h,
                       const IndexView &ghosted_local_2_global)
            : dims_(dims),
              local_elements_begin_(local_elements_begin),
              local_elements_end_(local_elements_end),
              box_(box),
              h_(h),
              ghosted_local_2_global_(ghosted_local_2_global) {}

        UTOPIA_INLINE_FUNCTION SizeType n_elements() const { return prod_dim(dims_, 0); }

        UTOPIA_INLINE_FUNCTION SizeType n_nodes() const { return prod_dim(dims_, 1); }

        UTOPIA_INLINE_FUNCTION SizeType n_local_elements() const {
            SizeType ret = local_elements_end_[0] - local_elements_begin_[0];
            for (int d = 1; d < Dim; ++d) {
                ret *= local_elements_end_[d] - local_elements_begin_[d];
            }

            return ret;
        }

        UTOPIA_INLINE_FUNCTION SizeType local_element_index_begin() const {
            return TensorMeshImpl<TensorMeshView>::local_element_index_begin(dims_, local_elements_begin_);
        }

        UTOPIA_INLINE_FUNCTION SizeType local_element_index_end() const {
            return TensorMeshImpl<TensorMeshView>::local_element_index_end(dims_, local_elements_end_);
        }

        template <class Array>
        UTOPIA_INLINE_FUNCTION SizeType element_tensor_to_linear_index(const Array &tensor_index) const {
            return TensorMeshImpl<TensorMeshView>::element_tensor_to_linear_index(dims_, tensor_index);
        }

        template <class Array>
        UTOPIA_INLINE_FUNCTION void element_linear_to_tensor_index(const SizeType &element_idx,
                                                                   Array &tensor_index) const {
            TensorMeshImpl<TensorMeshView>::element_linear_to_tensor_index(dims_, element_idx, tensor_index);
        }

        UTOPIA_INLINE_FUNCTION SizeType local_node_index_begin() const {
            return TensorMeshImpl<TensorMeshView>::local_node_index_begin(dims_, local_elements_begin_);
        }

        UTOPIA_INLINE_FUNCTION SizeType local_node_index_end() const {
            return TensorMeshImpl<TensorMeshView>::local_node_index_end(dims_, local_elements_end_);
        }

        template <class Array>
        UTOPIA_INLINE_FUNCTION void node_linear_to_tensor_index(const SizeType &node_idx, Array &tensor_index) const {
            return TensorMeshImpl<TensorMeshView>::node_linear_to_tensor_index(dims_, node_idx, tensor_index);
        }

        template <class Array>
        UTOPIA_INLINE_FUNCTION SizeType node_tensor_to_linear_index(const Array &tensor_index) const {
            return TensorMeshImpl<TensorMeshView>::node_tensor_to_linear_index(dims_, tensor_index);
        }

        template <class Array>
        UTOPIA_INLINE_FUNCTION static SizeType prod_dim(Array &arr, const SizeType offset = 0) {
            SizeType ret = arr[0] + offset;
            for (int d = 1; d < Dim; ++d) {
                ret *= arr[d] + offset;
            }

            return ret;
        }

        template <class Array>
        UTOPIA_INLINE_FUNCTION bool is_boundary_node(const Array &tensor_index) const {
            for (int i = 0; i < Dim; ++i) {
                if (tensor_index[i] == 0 || tensor_index[i] == dims_[i]) {
                    return true;
                }
            }

            return false;
        }

        // dir is element of
        UTOPIA_INLINE_FUNCTION bool is_node_ghost_layer(const int &dim, const int &dir) const {
            if (dir < 0) {
                SizeType idx = local_elements_begin_[dim];
                return idx != 0;
            } else {
                SizeType idx = local_elements_end_[dim];
                return idx != dims_[dim];
            }

            return false;
        }

        UTOPIA_INLINE_FUNCTION SizeType n_nodes_in_slice(const int &dim) const {
            SizeType ret = 1;
            for (int i = 0; i < Dim; ++i) {
                if (i == dim) continue;
                ret *= (local_elements_end_[i] - local_elements_end_[i]);
            }

            return ret;
        }

        const IndexView &dims() const { return dims_; }

    private:
        IndexView dims_;
        IndexView local_elements_begin_;
        IndexView local_elements_end_;
        Box<ScalarView> box_;
        ScalarView h_;
        IndexView ghosted_local_2_global_;
    };

    template <class Elem_, typename...>
    class Mesh {};

    template <class Elem_, class Comm, class ExecutionSpace, typename... Args>
    class Mesh<Elem_, Comm, ExecutionSpace, Uniform<Args...>> {
    public:
        using SizeType = std::size_t;
        using Elem = Elem_;

        static const int Dim = Elem::Dim;
        static const std::size_t UDim = Dim;
        static const int NNodesXElem = utopia::Power<Dim, 2>::value;

        static const int LEFT = -1;
        static const int RIGHT = 1;

        using Scalar = typename Elem::Scalar;
        using DualIndexView = Kokkos::DualView<SizeType *, ExecutionSpace>;
        using DualScalarView = Kokkos::DualView<Scalar *, ExecutionSpace>;
        using Dev = utopia::TpetraTraits::Device;
        using DualBoxView = utopia::Box<DualScalarView>;

        using ViewDevice = utopia::TensorMeshView<Elem, ExecutionSpace>;

        Mesh(const Comm &comm,
             const DualIndexView &dims,
             const DualIndexView &local_elements_begin,
             const DualIndexView &local_elements_end,
             const DualBoxView &box,
             const DualIndexView &ghosted_local_2_global = DualIndexView())
            : comm_(comm),
              dims_(dims),
              local_elements_begin_(local_elements_begin),
              local_elements_end_(local_elements_end),
              box_(box),
              h_("h", Dim),
              ghosted_local_2_global_(ghosted_local_2_global) {
            init();
        }

        Mesh(const Comm &comm,
             const std::array<SizeType, UDim> &dims,
             const std::array<SizeType, UDim> &local_elements_begin,
             const std::array<SizeType, UDim> &local_elements_end,
             const Box<std::array<Scalar, UDim>> &box,
             const DualIndexView &ghosted_local_2_global = DualIndexView())
            : comm_(comm),
              dims_("dims", Dim),
              local_elements_begin_("local_elements_begin", Dim),
              local_elements_end_("local_elements_end", Dim),
              box_(box),
              h_("h", Dim),
              ghosted_local_2_global_(ghosted_local_2_global) {
            auto h_dim = dims_.view_host();
            auto h_local_elements_begin = local_elements_begin_.view_host();
            auto h_local_elements_end = local_elements_end_.view_host();

            for (int d = 0; d < Dim; ++d) {
                h_dim[d] = dims[d];
                h_local_elements_begin[d] = local_elements_begin[d];
                h_local_elements_end[d] = local_elements_end[d];
            }

            init();
        }

        void init() {
            auto max = box_.max.view_host();
            auto min = box_.min.view_host();
            auto h = h_.view_host();
            auto dims = dims_.view_host();

            for (int d = 0; d < Dim; ++d) {
                h[d] = (max[d] - min[d]) / dims[d];
            }
        }

        Range element_range() const { return Range(local_element_index_begin(), local_element_index_end()); }

        SizeType local_element_index_begin() const {
            return TensorMeshImpl<Mesh>::local_element_index_begin(dims_.view_host(),
                                                                   local_elements_begin_.view_host());
        }

        SizeType local_element_index_end() const {
            return TensorMeshImpl<Mesh>::local_element_index_end(dims_.view_host(), local_elements_end_.view_host());
        }

        Range local_node_range() const { return Range(local_node_index_begin(), local_node_index_end()); }

        SizeType local_node_index_begin() const {
            return TensorMeshImpl<Mesh>::local_node_index_begin(dims_.view_host(), local_elements_begin_.view_host());
        }

        SizeType local_node_index_end() const {
            return TensorMeshImpl<Mesh>::local_node_index_end(dims_.view_host(), local_elements_end_.view_host());
        }

        template <class Array>
        SizeType element_tensor_to_linear_index(const Array &tensor_index) const {
            return TensorMeshImpl<Mesh>::element_tensor_to_linear_index(dims_.view_host(), tensor_index.view_host());
        }

        template <class Fun>
        void each_element(Fun fun) {
            auto view = view_device();
            Dev::parallel_for(
                element_range(), UTOPIA_LAMBDA(const SizeType &e_index) {
                    typename ViewDevice::Elem e;
                    view.elem(e_index, e);
                    fun(e_index, e);
                });
        }

        ViewDevice view_device() const {
            return ViewDevice(dims_.view_device(),
                              local_elements_begin_.view_device(),
                              local_elements_end_.view_device(),
                              {box_.min.view_device(), box_.max.view_device()},
                              h_.view_device(),
                              ghosted_local_2_global_.view_device());
        }

        inline SizeType n_nodes() const {
            auto dims = dims_.view_host();
            SizeType ret = dims[0] + 1;
            for (int d = 1; d < Dim; ++d) {
                ret *= dims[d] + 1;
            }

            return ret;
        }

        inline const Comm &comm() const { return comm_; }

        inline void describe(std::ostream &os = std::cout) const {
            auto h_dims = dims_.view_host();
            auto h_local_elements_begin = local_elements_begin_.view_host();
            auto h_local_elements_end = local_elements_end_.view_host();
            // h_box = box_.view_host();

            os << "dims:\n";
            for (int i = 0; i < Dim; ++i) {
                os << h_dims[i] << std::endl;
            }

            os << "elem_begin:\n";
            for (int i = 0; i < Dim; ++i) {
                os << h_local_elements_begin[i] << std::endl;
            }

            os << "elem_end:\n";
            for (int i = 0; i < Dim; ++i) {
                os << h_local_elements_end[i] << std::endl;
            }
        }

    private:
        Comm comm_;
        DualIndexView dims_;
        DualIndexView local_elements_begin_;
        DualIndexView local_elements_end_;
        DualBoxView box_;
        DualScalarView h_;
        DualIndexView ghosted_local_2_global_;

        template <class Array>
        static SizeType prod_dim(Array &arr, const SizeType offset = 0) {
            SizeType ret = arr[0] + offset;
            for (int d = 1; d < Dim; ++d) {
                ret *= arr[d] + offset;
            }

            return ret;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MESH_VIEW_HPP
