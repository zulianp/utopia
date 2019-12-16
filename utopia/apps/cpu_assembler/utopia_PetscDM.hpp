#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP


#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_Box.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_ArrayView.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Quadrature.hpp"
#include "utopia_UniformQuad4.hpp"

#include <array>
#include <memory>
#include <functional>

namespace utopia {

    template<int Dim>
    class PetscDMElements;

    template<int Dim>
    class PetscDMNodes;

    template<int Dim_>
    class PetscElem {
    public:
        static const int Dim = Dim_;
        using SizeType = PetscInt;
        using Scalar   = PetscScalar;
        using Point    = utopia::StaticVector<Scalar, Dim>;
        using Grad     = utopia::StaticVector<Scalar, Dim>;
        using NodeIndex = utopia::ArrayView<const SizeType>;

        virtual ~PetscElem() {}
        PetscElem()
        {}

        inline void set(const NodeIndex &nodes)
        {
            nodes_ = nodes;
        }

        inline const NodeIndex &nodes() const
        {
            return nodes_;
        }


        inline NodeIndex &nodes()
        {
            return nodes_;
        }

        inline SizeType node_id(const SizeType k) const
        {
            return nodes_[k];
        }

        inline const SizeType &node(const SizeType &i) const
        {
            return nodes_[i];
        }

        inline SizeType n_nodes() const
        {
            return nodes_.size();
        }

        inline constexpr static SizeType dim() { return Dim; }
        virtual Scalar fun(const SizeType &i, const Point &p) const = 0;

        inline void idx(const SizeType &idx)
        {
            idx_ = idx;
        }

        inline SizeType idx() const
        {
            return idx_;
        }

    private:
        NodeIndex nodes_;
        SizeType idx_;
    };

    template<int Dim>
    class PetscNode {
    public:
        PetscNode(const PetscDMNodes<Dim> &nodes, const SizeType &idx) : nodes_(nodes), idx_(idx) {}
        SizeType idx() const { return idx_; }
        bool is_ghost() const;

    private:
        const PetscDMNodes<Dim> &nodes_;
        SizeType idx_;
    };

    template<int Dim>
    class PetscDM final {
    public:
        static const std::size_t UDim = Dim;
        using SizeType  = PetscInt;
        using Scalar    = PetscScalar;
        using Point     = utopia::StaticVector<Scalar, Dim>;
        using Elem      = utopia::PetscElem<Dim>;
        using Node      = utopia::PetscNode<Dim>;
        using Device    = utopia::Device<PETSC>;
        using NodeIndex = typename Elem::NodeIndex;

        class Impl;


        PetscDM(
            const PetscCommunicator     &comm,
            const std::array<SizeType, UDim> &dims,
            const std::array<Scalar, UDim>   &box_min,
            const std::array<Scalar, UDim>   &box_max
        );

        PetscDM();
        ~PetscDM();

        void cell_point(const SizeType &idx, Point &translation);
        void cell_size(const SizeType &idx, Point &cell_size);

        void elem(const SizeType &idx, Elem &e) const;
        void nodes(const SizeType &idx, NodeIndex &nodes) const;
        void nodes_local(const SizeType &idx, NodeIndex &nodes) const;

        void create_matrix(PetscMatrix &mat);
        void create_vector(PetscVector &vec);

        // void create_local_matrix(PetscMatrix &mat);
        void create_local_vector(PetscVector &vec);

        // void local_to_global(PetscMatrix &local, PetscMatrix &global);
        void local_to_global(PetscVector &local, PetscVector &global);

        void describe() const;

        void each_element(const std::function<void(const Elem &)> &f);
        void each_node(const std::function<void(const Node &)> &f);
        void each_node_with_ghosts(const std::function<void(const Node &)> &f);

        inline static constexpr SizeType dim() { return Dim; }
        Range local_node_range() const;
        SizeType n_local_nodes_with_ghosts() const;

        // void dims(SizeType *arr) const;
        void box(Scalar *min, Scalar *max) const;

        void local_node_ranges(SizeType *begin, SizeType *end) const;
        void local_element_ranges(SizeType *begin, SizeType *end) const;

        Range local_element_range() const;

        template<class Array>
        void local_element_ranges(Array &begin, Array &end) const
        {
            SizeType temp_begin[3], temp_end[3];
            local_element_ranges(temp_begin, temp_end);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                begin[i] = temp_begin[i];
                end[i] = temp_end[i];
            }
        }

        template<class Array>
        void local_node_ranges(Array &begin, Array &end) const
        {
            SizeType temp_begin[3], temp_end[3];
            local_node_ranges(temp_begin, temp_end);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                begin[i] = temp_begin[i];
                end[i] = temp_end[i];
            }
        }

        template<class Array>
        void dims(Array &arr) const
        {
            std::array<SizeType, UDim> temp;
            dims(temp);
            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                arr[i] = temp[i];
            }
        }

        void dims(std::array<SizeType, UDim> &arr) const;

        template<class Array>
        void box(Array &min, Array &max) const
        {
            Scalar temp_min[3], temp_max[3];
            box(temp_min, temp_max);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                min[i] = temp_min[i];
                max[i] = temp_max[i];
            }
        }

        inline Impl &impl()
        {
            return *impl_;
        }

        inline const Impl &impl() const
        {
            return *impl_;
        }

        bool is_ghost(const SizeType &global_node_idx) const;

    private:
        std::unique_ptr<Impl> impl_;
    };

    class PetscUniformQuad4 final : public PetscElem<2> {
    public:
        using Super    = utopia::PetscElem<2>;
        using SizeType = Super::SizeType;
        using Scalar   = Super::Scalar;
        using Point    = Super::Point;
        using Grad     = Super::Grad;
        using MemType  = Uniform<>;
        static const int Dim = 2;
        static const int NNodes = 4;

        inline Scalar fun(const SizeType &i, const Point &p) const
        {
            return impl_.fun(i, p);
        }

        inline void node(const SizeType &i, Point &p) const
        {
            return impl_.node(i, p);
        }

        inline void grad(const int i, const Point &p, Grad &g) const
        {
           impl_.grad(i, p, g);
        }

        inline constexpr static bool is_affine()
        {
            return UniformQuad4<Scalar>::is_affine();
        }

        inline constexpr static Scalar reference_measure()
        {
            return UniformQuad4<Scalar>::reference_measure();
        }

        inline constexpr static int n_nodes()
        {
            return UniformQuad4<Scalar>::n_nodes();
        }

        inline void set(
            const StaticVector2<Scalar> &translation,
            const StaticVector2<Scalar> &h)
        {
            impl_.set(translation, h);
        }

        void describe(std::ostream &os = std::cout) const
        {
            auto &t = impl_.translation();
            os << t(0) << " " << t(1) << "\n";
        }

    private:
        UniformQuad4<Scalar> impl_;
    };

    template<typename View, std::size_t N>
    class Accessor<std::array<VectorView<View>, N>> {
    public:
        using Matrix = std::array<VectorView<View>, N>;
        using T = typename Traits<View>::Scalar;

        static const T &get(const Matrix &t, const std::size_t &i, const std::size_t &j)
        {
            return t[i][j];
        }

        static void set(Matrix &t, const std::size_t &i, const std::size_t &j, const T &val)
        {
            t[i][j] = val;
        }
    };

    template<>
    class Quadrature<PetscUniformQuad4, 2, 2> {
    public:
        using Scalar   = PetscUniformQuad4::Scalar;
        using SizeType = PetscUniformQuad4::SizeType;
        using Point    = PetscUniformQuad4::Point;
        using ViewDevice = const Quadrature &;

        static const int Order   = 2;
        static const int Dim     = 2;
        static const int NPoints = 6;


        inline static constexpr int n_points()
        {
            return NPoints;
        }

        inline static constexpr int dim()
        {
            return Dim;
        }

        void init()
        {
            Quad4Quadrature<Scalar, Order, Dim, NPoints>::get(points_, weights_);
        }

        template<class Point>
        inline void point(const int qp_idx, Point &p) const
        {
            p[0] = points_[qp_idx][0];
            p[1] = points_[qp_idx][1];
        }

        inline const Scalar &weight(const int qp_idx) const
        {
            return weights_[qp_idx];
        }

        Quadrature()
        {
            init();
        }

        inline ViewDevice &view_device() const
        {
            return *this;
        }

    private:
        std::array<Point, NPoints> points_;
        std::array<Scalar, NPoints> weights_;
    };



    template<class Elem_>
    class FunctionSpace<PetscDM<Elem_::Dim>, 1, Elem_> {
    public:
        static const int Dim = Elem_::Dim;
        using Mesh = utopia::PetscDM<Dim>;
        using Elem = Elem_;
        using MemType = typename Elem::MemType;
        using Scalar = typename Mesh::Scalar;
        using SizeType = typename Mesh::SizeType;
        using Point = typename Mesh::Point;

        using ViewDevice = FunctionSpace;
        using Device = typename Mesh::Device;
        using DofIndex = typename Mesh::NodeIndex;

        FunctionSpace(const std::shared_ptr<Mesh> &mesh)
        : mesh_(mesh)
        {}

        FunctionSpace(Mesh &mesh)
        : mesh_(utopia::make_ref(mesh))
        {}

        template<class Fun>
        void each_element(Fun fun)
        {
            mesh_->each_element(fun);
        }

        void elem(const SizeType &idx, Elem &e) const;

        void dofs(const SizeType &idx, DofIndex &dofs) const
        {
            mesh_->nodes(idx, dofs);
        }

        ViewDevice &view_device()
        {
            return *this;
        }

        Range local_element_range() const
        {
            return mesh_->local_element_range();
        }

        void create_matrix(PetscMatrix &mat)
        {
            mesh_->create_matrix(mat);
        }

        void create_vector(PetscVector &vec)
        {
            mesh_->create_vector(vec);
        }

    private:
        std::shared_ptr<Mesh> mesh_;

        FunctionSpace(const FunctionSpace &other)
        : mesh_(other.mesh_)
        {}
    };

    template<class Space>
    class BoundaryCondition {};


    template<class Elem, int Components>
    class BoundaryCondition<FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>> {
    public:
        using FunctionSpace = utopia::FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>;
        using Point    = typename FunctionSpace::Point;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar   = typename FunctionSpace::Scalar;
        using DofIndex = typename FunctionSpace::DofIndex;
        static const int Dim = FunctionSpace::Dim;

        enum SideSet {
            LEFT = 0,
            RIGHT = 1,
            BOTTOM = 2,
            TOP = 3,
            FRONT = 4,
            BACK = 5
        };

        BoundaryCondition(
            const FunctionSpace &space,
            SideSet side_set,
            const std::function<Scalar(const Point &)> &fun,
            const int component = 0)
        : space_(space), side_set_(side_set), fun_(fun), component_(component)
        {}

        void apply(PetscMatrix &mat, PetscVector &vec) const
        {

        }

        template<class ElementMatrix, class ElementVector>
        void apply(
            const Elem &e,
            const DofIndex &ind,
            ElementMatrix &mat, ElementVector &vec)
        {

        }

    private:
        const FunctionSpace &space_;
        SideSet side_set_;
        std::function<Scalar(const Point &)> fun_;
        int component_;
    };

}

#endif //UTOPIA_PETSC_DM_HPP
