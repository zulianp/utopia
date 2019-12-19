#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_Box.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_ArrayView.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Quadrature.hpp"
#include "utopia_UniformQuad4.hpp"
#include "utopia_UniformHex8.hpp"
#include "utopia_Temp.hpp"
#include "utopia_Writable.hpp"
#include "utopia_petsc_FE.hpp"
#include "utopia_DeviceView.hpp"

#include <array>
#include <memory>
#include <functional>

namespace utopia {

    class SideSet {
    public:
        using BoundaryIdType = int;
        inline static constexpr BoundaryIdType left() { return 1; }
        inline static constexpr BoundaryIdType right() { return 2; }
        inline static constexpr BoundaryIdType bottom() { return 3; }
        inline static constexpr BoundaryIdType top() { return 4; }
        inline static constexpr BoundaryIdType front() { return 5; }
        inline static constexpr BoundaryIdType back() { return 6; }
        inline static constexpr BoundaryIdType invalid() { return 0; }
    };

    template<int Dim>
    class PetscDMElements;

    template<int Dim>
    class PetscDMNodes;

    template<class Space>
    class DirichletBoundaryCondition {};

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

        PetscCommunicator &comm();
        const PetscCommunicator &comm() const;

        void cell_point(const SizeType &idx, Point &translation);
        void cell_size(const SizeType &idx, Point &cell_size);

        bool is_local_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const;
        bool is_node_on_boundary(const SizeType &idx, SideSet::BoundaryIdType b_id) const;
        void node(const SizeType &idx, Point &node) const;
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

        SizeType n_nodes() const;

        bool is_ghost(const SizeType &global_node_idx) const;
        bool is_boundary(const SizeType &global_node_idx) const;
        SideSet::BoundaryIdType boundary_id(const SizeType &global_node_idx) const;

    private:
        std::unique_ptr<Impl> impl_;
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
        using Vector = utopia::PetscVector;
        using Matrix = utopia::PetscMatrix;
        using Comm = utopia::PetscCommunicator;
        using DirichletBC = utopia::DirichletBoundaryCondition<FunctionSpace>;

        bool write(const Path &path, const PetscVector &x) const;
        // {
        //     PetscViewer       viewer;

        //     PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
        //     PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK);
        //     DMView(mesh_->raw_type(), viewer);
        //     PetscObjectSetName((PetscObject)raw_type(x), "x");
        //     VecView(raw_type(x), viewer);
        //     PetscViewerDestroy(&viewer);
        //     return false;
        // }

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

        bool is_boundary_dof(const SizeType &idx) const
        {
            return mesh_->is_boundary(idx);
        }

        SideSet::BoundaryIdType boundary_id(const SizeType &idx) const
        {
            return mesh_->boundary_id(idx);
        }

        inline static constexpr SizeType component(const SizeType &)
        {
            return 0;
        }

        const ViewDevice &view_device() const
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

        static DeviceView<PetscMatrix, 2> assembly_view_device(PetscMatrix &mat)
        {
            return DeviceView<PetscMatrix, 2>(mat, utopia::GLOBAL_ADD);
        }

        static DeviceView<PetscVector, 1> assembly_view_device(PetscVector &vec)
        {
            return DeviceView<PetscVector, 1>(vec, utopia::GLOBAL_ADD);
        }

        static DeviceView<const PetscVector, 1> assembly_view_device(const PetscVector &vec)
        {
            return DeviceView<const PetscVector, 1>(vec);
        }

        template<class ElementMatrix, class MatView>
        static void add_matrix(
            const Elem &e,
            const ElementMatrix &el_mat,
            MatView &mat)
        {
            const SizeType n_dofs = e.nodes().size();
            const auto &dofs = e.nodes();

            for(SizeType i = 0; i < n_dofs; ++i) {
                for(SizeType j = 0; j < n_dofs; ++j) {
                    mat.atomic_add(dofs[i], dofs[j], el_mat(i, j));
                }
            }
        }

        template<class ElementVector, class VecView>
        static void add_vector(
            const Elem &e,
            const ElementVector &el_vec,
            VecView &vec)
        {
            const SizeType n_dofs = e.nodes().size();
            const auto &dofs = e.nodes();

            for(SizeType i = 0; i < n_dofs; ++i) {
                vec.atomic_add(dofs[i], el_vec(i));
            }
        }

        //FIXME does not work in parallel
        template<class VectorView, class Values>
        static void coefficients(
            const Elem &e,
            const VectorView &vec,
            Values &values)
        {
            const SizeType n_dofs = e.nodes().size();
            const auto &dofs = e.nodes();

            for(SizeType i = 0; i < n_dofs; ++i) {
                values[i] = vec.get(dofs[i]);
            }
        }

        const Mesh &mesh() const
        {
            return *mesh_;
        }

        PetscCommunicator &comm()
        {
            return mesh_->comm();
        }

        const PetscCommunicator &comm() const
        {
            return mesh_->comm();
        }

        inline SizeType n_dofs() const
        {
            return mesh_->n_nodes();
        }

        template<class... Args>
        void emplace_dirichlet_condition(Args && ...args)
        {
            dirichlet_bcs_.push_back(utopia::make_unique<DirichletBC>(*this, std::forward<Args>(args)...));
        }

        void apply_constraints(PetscMatrix &mat, PetscVector &vec) const
        {
            for(const auto &bc : dirichlet_bcs_) {
                bc->apply(mat, vec);
            }
        }

        void apply_constraints(PetscMatrix &mat) const
        {
            for(const auto &bc : dirichlet_bcs_) {
                bc->apply(mat);
            }
        }

        void apply_constraints(PetscVector &vec) const
        {
            for(const auto &bc : dirichlet_bcs_) {
                bc->apply(vec);
            }
        }

        ///shallow copy
        FunctionSpace(const FunctionSpace &other)
        : mesh_(other.mesh_)
        {}

    private:
        std::shared_ptr<Mesh> mesh_;
        std::vector<std::shared_ptr<DirichletBC>> dirichlet_bcs_;
    };

    template<class Elem, int Components>
    class DirichletBoundaryCondition<FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>> {
    public:
        using FunctionSpace = utopia::FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>;
        using Point    = typename FunctionSpace::Point;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar   = typename FunctionSpace::Scalar;
        using DofIndex = typename FunctionSpace::DofIndex;
        static const int Dim = FunctionSpace::Dim;

        DirichletBoundaryCondition(
            const FunctionSpace &space,
            SideSet::BoundaryIdType side_set,
            const std::function<Scalar(const Point &)> &fun,
            const int component = 0)
        : space_(space), side_set_(side_set), fun_(fun), component_(component)
        {}

        void apply(PetscMatrix &mat, PetscVector &vec) const
        {
            using IndexSet = Traits<PetscVector>::IndexSet;
            IndexSet ind;
            ind.reserve(vec.local_size());

            auto r = vec.range();

            Write<PetscVector> w(vec, utopia::AUTO);

            Point p;
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(is_constrained_dof(i - r.begin())) {
                    ind.push_back(i);
                    space_.mesh().node(i/Components, p);
                    vec.set(i, fun_(p));
                }
            }

            mat.set_zero_rows(ind, 1.0);
        }

        void apply(PetscMatrix &mat) const
        {
            using IndexSet = Traits<PetscVector>::IndexSet;
            IndexSet ind;
            ind.reserve(mat.local_size().get(0));

            auto r = mat.row_range();

            Point p;
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(is_constrained_dof(i - r.begin())) {
                    ind.push_back(i);
                    space_.mesh().node(i/Components, p);
                }
            }

            mat.set_zero_rows(ind, 1.0);
        }

        void apply(PetscVector &vec) const
        {
            auto r = vec.range();

            Write<PetscVector> w(vec, utopia::AUTO);

            Point p;
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(is_constrained_dof(i - r.begin())) {
                    space_.mesh().node(i/Components, p);
                    vec.set(i, fun_(p));
                }
            }
        }

        void set_boundary_id(PetscVector &vec) const
        {
            auto r = vec.range();

            Write<PetscVector> w(vec, utopia::AUTO);
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(is_constrained_dof(i - r.begin())) {
                    vec.set(i, side_set_);
                }
            }
        }

        template<class ElementMatrix, class ElementVector>
        void apply(
            const Elem &e,
            const DofIndex &ind,
            ElementMatrix &mat, ElementVector &vec)
        {
            assert(false && "FIXME: use local node index or cartesian for identifying boundary nodes");

            const SizeType n_dofs = ind.size();
            Point p;

            for(SizeType i = 0; i < n_dofs; ++i) {
                if(!is_constrained_dof(ind[i])) continue;

                for(SizeType j = 0; j < n_dofs; ++j) {
                    mat(i, j) = (i == j);
                }

                e.node(i/Components, p);
                vec[i] = fun_(p);
            }
        }

    private:
        const FunctionSpace &space_;
        SideSet::BoundaryIdType side_set_;
        std::function<Scalar(const Point &)> fun_;
        int component_;

        bool is_constrained_dof(const SizeType &idx) const
        {
            if(space_.component(idx) != component_) return false;
            return space_.mesh().is_local_node_on_boundary(idx/Components, side_set_);
        }
    };

}

#endif //UTOPIA_PETSC_DM_HPP
