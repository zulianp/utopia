#ifndef UTOPIA_PETSC_NEUMANN_BOUNDARY_CONDITIONS_HPP
#define UTOPIA_PETSC_NEUMANN_BOUNDARY_CONDITIONS_HPP

#include "utopia_petsc_dma_FunctionSpace.hpp"

namespace utopia {

    template<class Elem, int Components>
    class NeumannBoundaryCondition<FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>> {
    public:
        using Mesh    = utopia::PetscDM<Elem::Dim>;
        using NodeIndex = typename Mesh::NodeIndex;
        using FunctionSpace = utopia::FunctionSpace<Mesh, Components, Elem>;
        using Point    = typename FunctionSpace::Point;
        using SizeType = typename FunctionSpace::SizeType;
        using Scalar   = typename FunctionSpace::Scalar;
        using Subspace = typename FunctionSpace::template Subspace<1>;
        using Device   = typename Subspace::Device;
        using ElemView = typename Subspace::ViewDevice::Elem;
        using IndexSet = Traits<PetscVector>::IndexSet;

        static const int Dim = Subspace::Dim;
        static const int NFunctions = Subspace::NFunctions;

        NeumannBoundaryCondition(const FunctionSpace &space)
        : space_(space), side_set_(SideSet::invalid()), component_(0)
        {}

        NeumannBoundaryCondition(
            const FunctionSpace &space,
            SideSet::BoundaryIdType side_set,
            const std::function<Scalar(const Point &)> &fun,
            const int component = 0)
        : space_(space), side_set_(SideSet::invalid()), fun_(fun), component_(component)
        {
            // init_element_list();
        }


        void apply(PetscVector &v) const
        {
            using Side     = typename Elem::Side;
            using SideView = typename Side::ViewDevice;
            using SideQuadrature = utopia::Quadrature<SideElem, 2, 1>;

            auto subspace   = space_.subspace(component_);

            SideQuadrature q;
            auto points      = subspace.side_points(q);
            auto shape       = subspace.side_shape(q);

            auto space_view  = subspace.view_device();
            auto v_view      = subspace.assembly_view(v);
            auto points_view = points.view_device();
            auto shape_view  = shape.view_device();

            Device::parallel_for(
                // subspace.boundary_element_range(side_set_),
                subspace.local_element_range(),
                UTOPIA_LAMBDA(const SizeType &i)
            {
                if(!space_view.on_boundary(i)) {
                    return;
                }

                StaticVector<Scalar, NFunctions> vec; vec.set(0.0);
                SideElemView e;

                auto p   = points_view.make(e);
                auto fun = shape_view.make(e);

                for(SizeType s = 0; s < Elem::NSides; ++s) {
                    if(!space_view.on_boundary(i, s, side_set_)) {
                        continue;
                    }

                    space_view.elem(i, s, e);

                    //assemble
                    for(SizeType k = 0; k < SideQuadrature::NPoints; ++k) {
                        auto dx_k  = dx(k);
                        auto p_k   = p(k);

                        if(!selector_ || selector_(p_k)) {

                            auto fun_k = fun_(p_k);

                            for(SizeType j = 0; j < Side::NFunctions; ++j) {
                                vec(j) += fun_k * fun(j, k) * dx_k;
                            }

                        }
                    }
                }

                space_view.add_vector(e, vec, v_view);
            });

        }

    private:
        const FunctionSpace &space_;
        SideSet::BoundaryIdType side_set_;

        std::function<bool(const Point &)>   selector_;
        std::function<Scalar(const Point &)> fun_;
        int component_;

        // IndexSet indices_;

        // void init_element_list()
        // {

        // }
    };

}

#endif //UTOPIA_PETSC_NEUMANN_BOUNDARY_CONDITIONS_HPP