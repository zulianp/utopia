#ifndef UTOPIA_PETSC_DIRICHLET_BOUNDARY_CONDITIONS_HPP
#define UTOPIA_PETSC_DIRICHLET_BOUNDARY_CONDITIONS_HPP

#include "utopia_petsc_dma_FunctionSpace.hpp"

namespace utopia {

    template<class Elem, int Components>
    class DirichletBoundaryCondition<FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>> {
    public:
        using FunctionSpace = utopia::FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>;
        using Point    = typename FunctionSpace::Point;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar   = typename FunctionSpace::Scalar;
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

        template<class DofIndex, class ElementMatrix, class ElementVector>
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


#endif //UTOPIA_PETSC_DIRICHLET_BOUNDARY_CONDITIONS_HPP
