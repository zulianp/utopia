#ifndef UTOPIA_TRILINOS_STRUCTURED_GRID_HPP
#define UTOPIA_TRILINOS_STRUCTURED_GRID_HPP

#include "utopia_trilinos.hpp"
#include "utopia_trilinos_FECrsGraph.hpp"
#include "utopia_Views.hpp"


namespace utopia {

    class FE {};
    class FD {};
    class FV {};

    template<
        class Elem_,
        typename...>
    class Mesh {};

    template<class Mesh, int NComponents, typename ...>
    class FunctionSpace {};

    template<typename T, int Order, int Dim = T::Dim, typename ...>
    class Quadrature {};

    template<typename...>
    class Uniform {};

    template<typename...>
    class Varying {};

    class RefQuad4 {
    public:
        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(p[0])
        {
            const auto x = p[0];
            const auto y = p[1];

            switch(i) {
                case 0: { return (1 - x) * (1 - y); }
                case 1: { return x * (1 - y); }
                case 2: { return x * y; }
                case 3: { return (1 - x) * y; }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return 0.0;
                }
            }
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION static void grad(const int i, const Point &p, Grad &g)
        {
            const auto x = p[0];
            const auto y = p[1];

            switch(i)
            {
                case 0:
                {
                    g[0] = y - 1.;
                    g[1] = x - 1.;
                    return;
                }
                case 1:
                {
                    g[0] = 1 - y;
                    g[1] = -x;
                    return;
                }
                case 2:
                {
                    g[0] = y;
                    g[1] = x;
                    return;
                }
                case 3:
                {
                    g[0] = -y;
                    g[1] = (1 - x);
                    return;
                }
                default:
                {
                    g[0] = 0.0;
                    g[1] = 0.0;
                    return;
                }
            }
        }

        template<typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION static void fun(const Point &p, Values &values)
        {
            const auto x = p[0];
            const auto y = p[1];

            values[0] = (1 - x) * (1 - y);
            values[1] = x * (1 - y);
            values[2] = x * y;
            values[3] = (1 - x) * y;
        }

        template<typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION static void grad(const Point &p, Values &values)
        {
            const auto x = p[0];
            const auto y = p[1];

            values(0, 0) = y - 1.;
            values(0, 1) = x - 1.;

            values(1, 0) = 1 - y;
            values(1, 1) = -x;

            values(2, 0) = y;
            values(2, 1) = x;

            values(3, 0) = -y;
            values(3, 1) = (1 - x);
        }
    };

    template<typename Scalar_>
    class UniformQuad4 final {
    public:
        using Scalar = Scalar_;
        using MemType = Uniform<>;
        using DiscretizationType = FE;
        static const int Dim = 2;
        static const int NNodes = 4;

        using NodeIndexView = utopia::ArrayView<std::size_t, NNodes>;

        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefQuad4::fun(i, p))
        {
          return RefQuad4::fun(i, p);
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const
        {
            RefQuad4::grad(i, p, g);
            g[0] /= h_[0];
            g[1] /= h_[1];
        }

        UTOPIA_INLINE_FUNCTION constexpr static bool is_affine()
        {
            return true;
        }

        UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure()
        {
            return 1.0;
        }

        // template<typename Point, typename Values>
        // UTOPIA_INLINE_FUNCTION static void fun(const Point &p, Values &values)
        // {
        //     RefQuad4::fun(p, values);
        // }

        // template<typename Point, typename Values>
        // UTOPIA_INLINE_FUNCTION void grad(const Point &p, Values &values) const
        // {
        //     RefQuad4::grad(p, values);
        //     for(int i = 0; i < 4; ++i)
        //     {
        //         values(i, 0) /= h_[0];
        //         values(i, 1) /= h_[1];
        //     }
        // }

        template<typename RefPoint, typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const
        {
            out[0] = in[0] * h_[0] + translation_[0];
            out[1] = in[1] * h_[1] + translation_[1];
        }

        UTOPIA_INLINE_FUNCTION UniformQuad4(const Scalar &hx, const Scalar &hy)
        {
            h_[0] = hx;
            h_[1] = hy;
        }

        UTOPIA_INLINE_FUNCTION UniformQuad4()
        {
            h_[0] = 0.0;
            h_[1] = 0.0;
        }

        template<class H>
        UTOPIA_INLINE_FUNCTION void set(
            const StaticVector2<Scalar> &translation,
            const H &h)
        {
            translation_(0) = translation(0);
            translation_(1) = translation(1);

            h_[0] = h[0];
            h_[1] = h[1];
        }

        UTOPIA_INLINE_FUNCTION NodeIndexView &nodes()
        {
            return nodes_;
        }

        UTOPIA_INLINE_FUNCTION const NodeIndexView &nodes() const
        {
            return nodes_;
        }

        UTOPIA_INLINE_FUNCTION const std::size_t &node(const std::size_t &i) const
        {
            return nodes_[i];
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_nodes()
        {
            return NNodes;
        }

    private:
        Scalar h_[2];
        StaticVector2<Scalar> translation_;
        NodeIndexView nodes_;
    };

    template<typename Scalar_>
    class Quadrature<UniformQuad4<Scalar_>, 2, 2> {
    public:
        using Scalar = Scalar_;
        using PointView   = Kokkos::View<Scalar **>;
        using WeightView  = Kokkos::View<Scalar *>;

        using DualPointView   = Kokkos::DualView<Scalar **>;
        using DualWeightView  = Kokkos::DualView<Scalar *>;

        static const int Order   = 2;
        static const int Dim     = 2;
        static const int NPoints = 6;

        UTOPIA_INLINE_FUNCTION static constexpr int n_points()
        {
            return NPoints;
        }

        UTOPIA_INLINE_FUNCTION static constexpr int dim()
        {
            return Dim;
        }

        UTOPIA_INLINE_FUNCTION const PointView &points() const
        {
            return points_;
        }

        template<class Point>
        UTOPIA_INLINE_FUNCTION void point(const int qp_idx, Point &p) const
        {
            p[0] = points_(qp_idx, 0);
            p[1] = points_(qp_idx, 1);
        }

        UTOPIA_INLINE_FUNCTION const Scalar &weight(const int qp_idx) const
        {
            return weights_[qp_idx];
        }

        UTOPIA_INLINE_FUNCTION const WeightView &weights() const
        {
            return weights_;
        }

        void init()
        {
            DualPointView points("Quadrature::points", NPoints, Dim);
            DualWeightView weights("Quadrature::weights", NPoints);

            auto host_points  = points.view_host();
            auto host_weights = weights.view_host();

            host_points(0, 0) = 0.5;
            host_points(0, 1) = 0.5;
            host_points(1, 0) = 0.9304589153964795245728880523899,
            host_points(1, 1) = 0.5;
            host_points(2, 0) = 0.72780186391809642112479237299488;
            host_points(2, 1) = 0.074042673347699754349082179816666;
            host_points(3, 0) = 0.72780186391809642112479237299488;
            host_points(3, 1) = 0.92595732665230024565091782018333;
            host_points(4, 0) = 0.13418502421343273531598225407969;
            host_points(4, 1) = 0.18454360551162298687829339850317;
            host_points(5, 0) = 0.13418502421343273531598225407969;
            host_points(5, 1) = 0.81545639448837701312170660149683;

            host_weights(0) = 0.28571428571428571428571428571428;
            host_weights(1) = 0.10989010989010989010989010989011;
            host_weights(2) = 0.14151805175188302631601261486295;
            host_weights(3) = 0.14151805175188302631601261486295;
            host_weights(4) = 0.16067975044591917148618518733485;
            host_weights(5) = 0.16067975044591917148618518733485;

            points_ = points.view_device();
            weights_ = weights.view_device();
        }

        Quadrature()
        {
            init();
        }

    private:
        PointView points_;
        WeightView weights_;
    };

    template<class Mesh, class MemType = typename Mesh::Elem::MemType>
    class FEValues {};

    template< int X, int Exponent >
    struct StaticPow {
        enum{ value = X*StaticPow<X, Exponent-1>::value };
    };

    template< int X >
    struct StaticPow< X, 1 > {
        enum{ value = X };
    };

    template<int Dim>
    class CellToNodeIndex {};


    template<>
    class CellToNodeIndex<2> {
    public:
        using SizeType = std::size_t;

        template<class Dims, class TensorIndexView, class NodeIndexView>
        UTOPIA_INLINE_FUNCTION static void build(
            const Dims &dims,
            const TensorIndexView &cell_tensor_idx,
            NodeIndexView &node_idx)
        {
            const SizeType i = cell_tensor_idx[0];
            const SizeType j = cell_tensor_idx[1];
            const SizeType n_p = dims[1] + 1;

            //flipped for dof-consistency and ccw local

            node_idx[0] = i * n_p + j;
            node_idx[3] = i * n_p + (j + 1);
            node_idx[2] = (i + 1) * n_p + (j + 1);
            node_idx[1] = (i + 1) * n_p + j;
        }
    };

    template<class Scalar, int Dim>
    class Box {
    public:
        ArrayView<Scalar, Dim> min;
        ArrayView<Scalar, Dim> max;
    };

    template<class Elem_, class Comm, class ExecutionSpace, typename...Args>
    class Mesh<Elem_, Comm, ExecutionSpace, Uniform<Args...>> {
    public:
        using SizeType = std::size_t;
        using Elem     = Elem_;

        static const int Dim = Elem::Dim;
        static const int NNodesXElem = utopia::StaticPow<Dim, 2>::value;

        using Scalar          = typename Elem::Scalar;
        using NodeIndexView   = utopia::ArrayView<SizeType, NNodesXElem>;
        using TensorIndexView = utopia::ArrayView<SizeType, Dim>;
        using IndexView1      = Kokkos::View<SizeType *, ExecutionSpace>;
        using ScalarView1     = Kokkos::View<Scalar   *, ExecutionSpace>;
        using Dev             = utopia::TpetraTraits::Device;


        template<class Array>
        UTOPIA_INLINE_FUNCTION void nodes(
            const SizeType &element_idx,
            Array &node_idx) const
        {
            TensorIndexView tensor_index;
            element_linear_to_tensor_index(element_idx, tensor_index);
            CellToNodeIndex<Dim>::build(dims_, tensor_index, node_idx);
        }

        template<class Point>
        UTOPIA_INLINE_FUNCTION void point(const SizeType &point_idx, Point &point)
        {
            TensorIndexView tensor_index;
            node_linear_to_tensor_index(point_idx, tensor_index);

            for(int d = 0; d < Dim; ++d) {
                point[d] = tensor_index[d] * h_[d] + box_.min[d];
            }
        }

        UTOPIA_INLINE_FUNCTION void elem(const SizeType &element_idx, Elem &elem) const
        {
            TensorIndexView tensor_index;
            element_linear_to_tensor_index(element_idx, tensor_index);

            StaticVector2<Scalar> translation;
            for(int d = 0; d < Dim; ++d) {
                translation(d) = tensor_index[d] * h_[d] + box_.min[d];
            }

            elem.set(translation, h_);
            nodes(element_idx, elem.nodes());
        }

        //FIXME this
        Mesh(
            const Comm &comm,
            const IndexView1 &dims,
            const IndexView1 &local_elements_begin,
            const IndexView1 &local_elements_end,
            const Box<Scalar, Dim> &box
            ) : comm_(comm),
                dims_(dims),
                local_elements_begin_(local_elements_begin),
                local_elements_end_(local_elements_end),
                box_(box),
                h_("h", Dim)
        {
            for(int d = 0; d < Dim; ++d) {
                h_[d] = (box_.max[d] - box_.min[d])/dims_[d];
            }
        }

        UTOPIA_INLINE_FUNCTION SizeType n_elements() const
        {
           return prod_dim(dims_, 0);
        }

        UTOPIA_INLINE_FUNCTION SizeType n_local_elements() const
        {
            SizeType ret = local_elements_end_[0] - local_elements_begin_[0];
            for(int d = 1; d < Dim; ++d) {
                ret *= local_elements_end_[d] - local_elements_begin_[d];
            }

            return ret;
        }

        UTOPIA_INLINE_FUNCTION SizeType n_nodes() const
        {
            return prod_dim(dims_, 1);
        }

        class TensorIndexIterator {
        public:
            UTOPIA_INLINE_FUNCTION TensorIndexIterator &operator++()
            {
                for(int d = Dim-1; d > 0; --d) {
                    index[d] += 1;
                    if(index[d] < end[d]) break;

                    index[d] = 0;
                    index[d-1] += 1;
                }

                return *this;
            }

            UTOPIA_INLINE_FUNCTION operator bool() const
            {
                for(int d = 0; d < Dim; ++d) {
                    if(index[d] < end[d]) return true;
                }

                return false;
            }

            const TensorIndexView & operator *() const
            {
                return index;
            }

            TensorIndexView index;
            TensorIndexView end;
        };

        UTOPIA_INLINE_FUNCTION SizeType local_elements_index_begin() const
        {
            return element_tensor_to_linear_index(local_elements_begin_);
        }

        UTOPIA_INLINE_FUNCTION SizeType local_elements_index_end() const
        {
            SizeType result = local_elements_end_[0] -1;

            for(int i = 1; i < Dim; ++i){
                result *= dims_[i];
                result += local_elements_end_[i] -1;
            }

            return result + 1;
        }

        Range local_element_range() const
        {
            return Range(local_elements_index_begin(), local_elements_index_end());
        }

        template<class Fun>
        void each_element(Fun fun)
        {
            Dev::parallel_for(local_element_range(), UTOPIA_LAMBDA(const SizeType &e_index) {
                Elem e;
                elem(e_index, e);
                fun(e_index, e);
            });
        }

    private:
        Comm comm_;
        IndexView1  dims_;
        IndexView1  local_elements_begin_;
        IndexView1  local_elements_end_;

        //FIXME
        Box<Scalar, Dim> box_;
        ScalarView1 h_;


        template<class Array>
        UTOPIA_INLINE_FUNCTION void element_linear_to_tensor_index(
            const SizeType &element_idx,
            Array &tensor_index
            ) const
        {
            SizeType current = element_idx;
            const int last = Dim - 1;

            for(int i = last; i >= 0; --i) {
                const int next = current / dims_[i];
                tensor_index[i] = current - next * dims_[i];
                current = next;
            }

            UTOPIA_DEVICE_ASSERT(element_tensor_to_linear_index(tensor_index) == element_idx);
        }

        template<class Array>
        UTOPIA_INLINE_FUNCTION SizeType element_tensor_to_linear_index(const Array &tensor_index) const
        {
            SizeType result = tensor_index[0];

            for(int i = 1; i < Dim; ++i){
                result *= dims_[i];
                result += tensor_index[i];
            }

            return result;
        }

        template<class Array>
        UTOPIA_INLINE_FUNCTION void node_linear_to_tensor_index(
            const SizeType &node_idx,
            Array &tensor_index
            ) const
        {
            SizeType current = node_idx;
            const int last = Dim - 1;

            for(int i = last; i >= 0; --i) {
                const SizeType next = current / (dims_[i] + 1);
                tensor_index[i] = current - next * (dims_[i] + 1);
                current = next;
            }

            UTOPIA_DEVICE_ASSERT(node_tensor_to_linear_index(tensor_index) == node_idx);
        }

        template<class Array>
        UTOPIA_INLINE_FUNCTION SizeType node_tensor_to_linear_index(const Array &tensor_index) const
        {
            SizeType result = tensor_index[0];

            for(SizeType i = 1; i < Dim; ++i){
                result *= dims_[i] + 1;
                result += tensor_index[i];
            }

            return result;
        }

        template<class Array>
        UTOPIA_INLINE_FUNCTION static SizeType prod_dim(Array &arr, const SizeType offset = 0)
        {
            SizeType ret = arr[0] + offset;
            for(int d = 1; d < Dim; ++d) {
                ret *= arr[d] + offset;
            }

            return ret;
        }
    };


    template<class Mesh_, int NComponents>
    class FunctionSpace<Mesh_, NComponents> {
    public:
        using Mesh = Mesh_;
        Mesh mesh_;
    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class Jacobian {};

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PhysicalPoint {
    public:
        static const int Dim = Elem::Dim;

        UTOPIA_INLINE_FUNCTION PhysicalPoint(const Elem &elem, const Quadrature &q)
        : elem_(elem), q_(q)
        {}

        template<class Point>
        UTOPIA_INLINE_FUNCTION void get(const int qp_idx, Point &p) const
        {
            Point temp;
            q_.point(qp_idx, temp);
            elem_.point(temp, p);
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return q_.n_points();
        }

    private:
        const Elem &elem_;
        const Quadrature &q_;
    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PhysicalGradient {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point  = utopia::StaticVector<Scalar, Dim>;

        UTOPIA_INLINE_FUNCTION PhysicalGradient(const Elem &elem, const Quadrature &q)
        : elem_(elem), q_(q)
        {}

        template<class Grad>
        UTOPIA_INLINE_FUNCTION void get(const int fun_num, const int qp_idx, Grad &g) const
        {
            Point temp;
            q_.point(qp_idx, temp);
            elem_.grad(fun_num, temp, g);
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return q_.n_points();
        }

    private:
        const Elem &elem_;
        const Quadrature &q_;
    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class ShapeFunction {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point  = utopia::StaticVector<Scalar, Dim>;

        UTOPIA_INLINE_FUNCTION ShapeFunction(const Elem &elem, const Quadrature &q)
        : elem_(elem), q_(q)
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const int fun_num, const int qp_idx) const
        {
            return get(fun_num, qp_idx);
        }

        UTOPIA_INLINE_FUNCTION Scalar get(const int fun_num, const int qp_idx) const
        {
            Point temp;
            q_.point(qp_idx, temp);
            return elem_.fun(fun_num, temp);
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return q_.n_points();
        }

    private:
        const Elem &elem_;
        const Quadrature &q_;
    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class Differential {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
        using Point  = utopia::StaticVector<Scalar, Dim>;

        UTOPIA_INLINE_FUNCTION Differential(const Elem &elem, const Quadrature &q)
        : elem_(elem), q_(q)
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const int qp_idx) const
        {
            return get(qp_idx);
        }

        UTOPIA_INLINE_FUNCTION Scalar get(const int qp_idx) const
        {
            if(elem_.is_affine()) {
                return q_.weight(qp_idx) * elem_.reference_measure();
            } else {
                UTOPIA_DEVICE_ASSERT(false);
                return -1.0;
            }
        }

        UTOPIA_INLINE_FUNCTION std::size_t size() const
        {
            return q_.n_points();
        }

    private:
        const Elem &elem_;
        const Quadrature &q_;
    };


    template<class Elem_, class Comm, class ExecutionSpace_, typename...Args>
    class FunctionSpace< Mesh<Elem_, Comm, ExecutionSpace_, Uniform<Args...>>, 1> {
    public:
        using Mesh = utopia::Mesh<Elem_, Comm, ExecutionSpace_, Uniform<Args...>>;
        using ExecutionSpace = ExecutionSpace_;
        using SizeType = typename Mesh::SizeType;
        using DofIndex = utopia::ArrayView<std::size_t, Elem_::NNodes>;

        template<class Array>
        UTOPIA_INLINE_FUNCTION void dofs(const SizeType &element_idx, Array &indices) const
        {
            mesh_.nodes(element_idx, indices);
        }

        template<class Fun>
        void each_element(Fun fun)
        {
            // Dev::parallel_for(local_element_range(), UTOPIA_LAMBDA(const SizeType &e_index) {
            //     Elem e;
            //     elem(e_index, e);
            //     fun(e_index, e);
            // });

            mesh_.each_element(fun);
        }

        FunctionSpace(const Mesh &mesh) : mesh_(mesh) {}
    private:
        Mesh mesh_;
    };

    // template<class Elem, typename...Args>
    // class FEValues<Mesh<Elem>, Uniform<Args...>> {
    // public:
    //     using Scalar   = typename Elem::Scalar;
    //     using Mesh     = utopia::Mesh<Elem>;
    //     using FunView  = Kokkos::DualView<Scalar **>;
    //     using GradView = Kokkos::DualView<Scalar ***>;

    //     template<typename...GridArgs>
    //     void init(const FunctionSpace<Mesh> &space)
    //     {

    //     }

    //     // void
    // };

    using Grid2d = utopia::Mesh<UniformQuad4<double>, TrilinosCommunicator, TpetraVector::vector_type::execution_space, Uniform<>>;
}

#endif
