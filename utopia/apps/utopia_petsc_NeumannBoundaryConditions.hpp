#ifndef UTOPIA_PETSC_NEUMANN_BOUNDARY_CONDITIONS_HPP
#define UTOPIA_PETSC_NEUMANN_BOUNDARY_CONDITIONS_HPP

#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_std_function.hpp"

namespace utopia {

    template<class Space, typename...>
    class NeumannBoundaryCondition {};


    //!! Only for petsc
    template<class Elem, int Components>
    class NeumannBoundaryCondition<FunctionSpace<PetscDM<Elem::Dim>, Components, Elem>> : public Configurable {
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
        using SideView = typename ElemView::Side;
        using IndexSet = Traits<PetscVector>::IndexSet;

        static const int Dim = Subspace::Dim;
        static const int NFunctions = ElemView::NFunctions;

        NeumannBoundaryCondition(const FunctionSpace &space)
        : space_(space), side_set_(SideSet::invalid()), component_(0)
        {}

        NeumannBoundaryCondition(
            const FunctionSpace &space,
            const utopia::function<bool(const Point &)>   &selector,
            const utopia::function<Scalar(const Point &)> &fun,
            const int component = 0)
        : space_(space), side_set_(SideSet::invalid()), selector_(selector), fun_(fun), component_(component)
        {}

        static utopia::function<bool(const Point &)> make_boundary_selector(
            const SideSet::BoundaryIdType &side_set,
            const Point &box_min,
            const Point &box_max
        )
        {
            switch(side_set) {
                case SideSet::left():
                {
                    const Scalar min_x = box_min[0];
                    return UTOPIA_LAMBDA(const Point &x) -> bool {
                        return device::approxeq(min_x, x[0], device::epsilon<Scalar>());
                    };

                    break;
                }
                case SideSet::right():
                {
                    const Scalar max_x = box_max[0];
                    return UTOPIA_LAMBDA(const Point &x) -> bool {
                        return device::approxeq(max_x, x[0], device::epsilon<Scalar>());
                    };

                    break;
                }
                case SideSet::bottom():
                {
                    const Scalar min_y = box_min[1];
                    return UTOPIA_LAMBDA(const Point &x) -> bool {
                        return device::approxeq(min_y, x[1], device::epsilon<Scalar>());
                    };

                    break;
                }
                case SideSet::top():
                {
                    const Scalar max_y = box_max[1];
                    return UTOPIA_LAMBDA(const Point &x) -> bool {
                        return device::approxeq(max_y, x[1], device::epsilon<Scalar>());
                    };

                    break;
                }
                case SideSet::back():
                {
                    const Scalar min_z = box_min[2];
                    return UTOPIA_LAMBDA(const Point &x) -> bool {
                        return device::approxeq(min_z, x[2], device::epsilon<Scalar>());
                    };

                    break;
                }
                case SideSet::front():
                {
                    const Scalar max_z = box_max[2];
                    return UTOPIA_LAMBDA(const Point &x) -> bool {
                        return device::approxeq(max_z, x[2], device::epsilon<Scalar>());
                    };

                    break;
                }
                default:
                {
                    assert(false);
                    return nullptr;
                    break;
                }
            }
        }

        NeumannBoundaryCondition(
            const FunctionSpace &space,
            SideSet::BoundaryIdType side_set,
            const utopia::function<Scalar(const Point &)> &fun,
            const int component = 0)
        : space_(space), side_set_(side_set), fun_(fun), component_(component)
        {
            auto &&box_min = space.mesh().box_min();
            auto &&box_max = space.mesh().box_max();
            selector_ = make_boundary_selector(side_set_, box_min, box_max);
        }

        void read(Input &in) override
        {
            std::string side_name = "";
            Scalar value = 0.0;

            in.get("side", side_name);
            in.get("value", value);
            in.get("component",  component_);

            side_set_ = SideSet::from_name(side_name);

            if(side_set_ != SideSet::invalid()) {
                auto &&box_min = space_.mesh().box_min();
                auto &&box_max = space_.mesh().box_max();
                selector_ = make_boundary_selector(side_set_, box_min, box_max);

                fun_ = UTOPIA_LAMBDA(const Point &) -> Scalar {
                    return value;
                };

            } else {
                assert(false);
            }
        }

        void apply(PetscVector &v) const
        {
            using Side     = typename Elem::Side;
            using SideQuadrature = utopia::Quadrature<Side, 2, Dim-1>;

            auto subspace   = space_.subspace(component_);

            SideQuadrature q;

            // q.describe();
            auto space_view  = subspace.view_device();
            auto v_view      = subspace.assembly_view_device(v);
            auto points_view = subspace.side_points_device(q);
            auto shape_view  = subspace.side_shape_device(q);
            auto dx_view     = subspace.side_differential_device(q);

            // subspace.each_boundary_element(side_set_, ...
            Device::parallel_for(
                subspace.element_range(),
                UTOPIA_LAMBDA(const SizeType &i)
            {
                if(!space_view.on_boundary(i)) {
                    return;
                }

                ArrayView<SizeType, Side::NFunctions> idx;
                StaticVector<Scalar, NFunctions> vec; vec.set(0.0);
                ElemView vol_e;
                SideView e;
                Point p_k;

                space_view.elem(i, vol_e);

                bool assembled = false;
                for(SizeType s = 0; s < Elem::NSides; ++s) {
                //TODO
                //     if(!space_view.on_boundary(vol_e, s, side_set_)) {
                //         continue;
                //     }

                    vol_e.side(s, e);

                    e.centroid(p_k);

                    if(!selector_(p_k)) continue;

                    assembled = true;

                    vol_e.side_idx(s, idx);

                    auto p   = points_view.make(e);
                    auto fun = shape_view.make(e);
                    auto dx  = dx_view.make(e);

                    for(SizeType k = 0; k < SideQuadrature::NPoints; ++k) {
                        auto dx_k = dx(k);
                        p.get(k, p_k);

                        auto fun_k = fun_(p_k);
                        for(SizeType j = 0; j < SideView::NFunctions; ++j) {
                            vec(idx[j]) += fun_k * fun(j, k) * dx_k;
                        }
                    }
                }

                if(assembled) {
                    space_view.add_vector(vol_e, vec, v_view);
                }
            });
        }


    private:
        const FunctionSpace &space_;
        SideSet::BoundaryIdType side_set_;

        utopia::function<bool(const Point &)>   selector_;
        utopia::function<Scalar(const Point &)> fun_;

        int component_;
    };

}

#endif //UTOPIA_PETSC_NEUMANN_BOUNDARY_CONDITIONS_HPP