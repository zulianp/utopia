#ifndef UTOPIA_PHASE_FIELD_HPP
#define UTOPIA_PHASE_FIELD_HPP

#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"

namespace utopia {

    template<class FunctionSpace>
    static void compute_strain_energy_splitting(
        const FunctionSpace &space,
        const typename FunctionSpace::Vector &displacement)
    {
        static const int Dim = FunctionSpace::Dim;
        static const int NVars = Dim;

        using Elem           = typename FunctionSpace::Elem;
        using ElemView       = typename FunctionSpace::ViewDevice::Elem;
        using Mesh           = typename FunctionSpace::Mesh;
        using SizeType       = typename Mesh::SizeType;
        using Scalar         = typename Mesh::Scalar;
        using Quadrature     = utopia::Quadrature<Elem, 2>;
        using Dev            = typename FunctionSpace::Device;
        using VectorD        = utopia::StaticVector<Scalar, Dim>;

        Quadrature q;
        PrincipalStrains<FunctionSpace, Quadrature> strain(space, q);
        strain.update(displacement);

        Differential<FunctionSpace, Quadrature> differential(space, q);

        Scalar mu = 1.0, lambda = 1.0;

        auto split_p = UTOPIA_LAMBDA(const Scalar &x) {
            return (device::abs(x) + x)/2;
        };

        auto split_m = UTOPIA_LAMBDA(const Scalar &x) {
            return (device::abs(x) - x)/2;
        };

        //end of host code

        StaticVector<Scalar, 2> energy;
        energy.set(0.0);
        {
            //beginning of device code
            auto strain_view = strain.view_device();
            auto space_view  = space.view_device();
            auto differential_view = differential.view_device();

            Dev::parallel_reduce(
                space.local_element_range(),
                UTOPIA_LAMBDA(const SizeType &i)
            {
                VectorD v;
                StaticVector<Scalar, 2> e_energy;

                ElemView e;
                space_view.elem(i, e);

                auto el_strain = strain_view.make(e);
                auto dx        = differential_view.make(e);

                const SizeType n_qp = el_strain.values.size();

                Scalar energy_p = 0;
                Scalar energy_m = 0;

                for(SizeType qp = 0; qp < n_qp; ++qp) {

                    Scalar trace_e_p = 0.0;
                    Scalar trace_e_n = 0.0;

                    Scalar inner_p = 0.0;
                    Scalar inner_n = 0.0;

                    //compute splitted quantities (inverting the sums)
                    for(int d = 0; d < Dim; ++d) {
                        auto e_val = el_strain.values[qp][d];

                        el_strain.vectors[qp].col(d, v);
                        auto outer_v = outer(v, v);
                        auto trace_v = trace(outer_v);
                        auto inner_v = inner(outer_v, outer_v);

                        auto eig_p = split_p(e_val);
                        auto eig_n = split_m(e_val);

                        trace_e_p += eig_p * trace_v;
                        trace_e_n += eig_n * trace_v;

                        inner_p += inner_v * eig_p * eig_p;
                        inner_n += inner_v * eig_n * eig_n;
                    }

                    energy_p += ((0.5 * lambda * trace_e_p * trace_e_p) +
                               (mu * inner_p)) * dx(qp);


                    energy_m += ((0.5 * lambda * trace_e_n * trace_e_n) +
                               (mu * inner_n)) * dx(qp);
                }

                e_energy[0] = energy_m;
                e_energy[1] = energy_p;
                return e_energy;
            }, energy);
        }

        space.comm().sum(2, &energy[0]);

        if(space.comm().rank() == 0) {
            disp(energy);
        }
    }
}
#endif