#ifndef UTOPIA_PETSC_DIRICHLET_BOUNDARY_CONDITIONS_HPP
#define UTOPIA_PETSC_DIRICHLET_BOUNDARY_CONDITIONS_HPP

#include "utopia_DeviceView.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_SideSets.hpp"

// petsc
#include "utopia_petsc.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"

namespace utopia {

    template <class FunctionSpace>
    class DirichletBoundaryCondition {
        // };

        // template<class Elem, int NComponents>
        // class DirichletBoundaryCondition<FunctionSpace<PetscDM<Elem::Dim>, NComponents, Elem>> {
    public:
        // using Mesh    = utopia::PetscDM<Elem::Dim>;

        using Mesh = typename FunctionSpace::Mesh;
        using NodeIndex = typename Mesh::NodeIndex;
        // using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Elem>;

        using Elem = typename FunctionSpace::Shape;
        using Point = typename FunctionSpace::Point;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar = typename FunctionSpace::Scalar;
        using Subspace = typename FunctionSpace::template Subspace<1>;
        using Device = typename Subspace::Device;
        using ElemView = typename Subspace::ViewDevice::Elem;
        using IndexSet = Traits<PetscVector>::IndexSet;

        static const int Dim = Subspace::Dim;
        static const int NComponents = FunctionSpace::NComponents;

        DirichletBoundaryCondition(const FunctionSpace &space)
            : space_(space), side_set_(SideSet::invalid()), component_(0) {}

        DirichletBoundaryCondition(const FunctionSpace &space,
                                   SideSet::BoundaryIdType side_set,
                                   const std::function<Scalar(const Point &)> &fun,
                                   const int component = 0)
            : space_(space), side_set_(side_set), fun_(fun), component_(component) {
            init_constraints_marker();
        }

        void apply(PetscMatrix &mat, PetscVector &vec) const {
            apply(mat);
            apply(vec);
        }

        void apply(PetscMatrix &mat) const { mat.set_zero_rows(indices_, 1.0); }

        void apply(PetscVector &v) const {
            // FIXME should just use the dofmap not the nodes

            auto r = v.range();

            auto subspace = space_.subspace(component_);

            auto space_view = subspace.view_device();
            auto v_view = utopia::view_device(v);

            Device::parallel_for(
                subspace.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                    ElemView e;
                    space_view.elem(i, e);

                    NodeIndex nodes;
                    space_view.mesh().nodes(i, nodes);

                    const SizeType n_nodes = nodes.size();

                    const SizeType nc = subspace.mesh().n_components();

                    Point p;
                    for (SizeType i = 0; i < n_nodes; ++i) {
                        const SizeType node_id = nodes[i];
                        auto idx = node_id * nc + component_;
                        e.node(i, p);

                        // FIXME
                        if (r.inside(idx) && is_constrained_dof(idx - r.begin())) {
                            v_view.set(idx, fun_(p));
                        }
                    }
                });
        }

        void apply_zero(PetscVector &vec) const { apply_val(vec, 0.0); }

        void apply_val(PetscVector &vec, const Scalar val) const {
            Write<PetscVector> w(vec, utopia::AUTO);
            for (auto i : indices_) {
                vec.set(i, val);
            }
        }

        void copy(const PetscVector &in, PetscVector &vec) const {
            Write<PetscVector> w(vec, utopia::AUTO);
            Read<PetscVector> r_in(in);

            for (auto i : indices_) {
                const Scalar val = in.get(i);
                vec.set(i, val);
            }
        }

        void set_boundary_id(PetscVector &vec) const {
            Write<PetscVector> w(vec, utopia::AUTO);

            for (auto i : indices_) {
                vec.set(i, side_set_);
            }
        }

        template <class DofIndex, class ElementMatrix, class ElementVector>
        void apply(const Elem &e, const DofIndex &ind, ElementMatrix &mat, ElementVector &vec) {
            assert(false && "FIXME: use local node index or cartesian for identifying boundary nodes");

            const SizeType n_dofs = ind.size();
            Point p;

            for (SizeType i = 0; i < n_dofs; ++i) {
                if (!is_constrained_dof(ind[i])) continue;

                for (SizeType j = 0; j < n_dofs; ++j) {
                    mat(i, j) = (i == j);
                }

                e.node(i / NComponents, p);
                vec[i] = fun_(p);
            }
        }

        void init_from(const DirichletBoundaryCondition &other) {
            side_set_ = other.side_set_;
            fun_ = other.fun_;
            component_ = other.component_;
            init_constraints_marker();
        }

    private:
        const FunctionSpace &space_;
        SideSet::BoundaryIdType side_set_;
        std::function<Scalar(const Point &)> fun_;
        int component_;

        IndexSet indices_;

        bool is_constrained_dof(const SizeType &idx) const {
            if (space_.component(idx) != component_) return false;
            return space_.mesh().is_node_on_boundary(idx / NComponents, side_set_);
        }

        void init_constraints_marker() {
            auto r = space_.dof_range();

            SizeType n_constrained_dofs = 0;
            for (auto i = r.begin(); i < r.end(); ++i) {
                if (is_constrained_dof(i - r.begin())) {
                    ++n_constrained_dofs;
                }
            }

            indices_.reserve(n_constrained_dofs);

            for (auto i = r.begin(); i < r.end(); ++i) {
                if (is_constrained_dof(i - r.begin())) {
                    indices_.push_back(i);
                }
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_DIRICHLET_BOUNDARY_CONDITIONS_HPP
