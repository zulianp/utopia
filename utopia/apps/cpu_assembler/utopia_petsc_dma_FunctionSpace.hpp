#ifndef UTOPIA_PETSC_DMA_FUNCTIONSPACE_HPP
#define UTOPIA_PETSC_DMA_FUNCTIONSPACE_HPP

#include "utopia_PetscDM.hpp"

namespace utopia {

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

        void dofs_local(const SizeType &idx, DofIndex &dofs) const
        {
            assert(idx < mesh_->n_elements());
            mesh_->nodes_local(idx, dofs);
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

        template<class VectorView, class Values>
        void local_coefficients(
            const Elem &e,
            const VectorView &vec,
            Values &values) const
        {
            DofIndex dofs;
            dofs_local(e.idx(), dofs);
            const SizeType n = dofs.size();

            for(SizeType i = 0; i < n; ++i) {
                assert(dofs[i] < n_dofs());
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

        void apply_zero_constraints(PetscVector &vec) const
        {
            for(const auto &bc : dirichlet_bcs_) {
                bc->apply_zero(vec);
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

        void apply_zero(PetscVector &vec) const
        {
            auto r = vec.range();

            Write<PetscVector> w(vec, utopia::AUTO);
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(is_constrained_dof(i - r.begin())) {
                    vec.set(i, 0.0);
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

#endif //UTOPIA_PETSC_DMA_FUNCTIONSPACE_HPP
