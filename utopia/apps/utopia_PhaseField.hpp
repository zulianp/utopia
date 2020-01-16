#ifndef UTOPIA_PHASE_FIELD_HPP
#define UTOPIA_PHASE_FIELD_HPP

#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_FEFunction.hpp"

namespace utopia {

    template<class FunctionSpace, int Dim = FunctionSpace::Dim>
    class PhaseFieldForBrittleFractures {
    public:
        using Scalar   = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector   = typename FunctionSpace::Vector;
        using Matrix   = typename FunctionSpace::Matrix;
        using Device   = typename FunctionSpace::Device;

        using USpace   = typename FunctionSpace::template Subspace<Dim>;
        using CSpace   = typename FunctionSpace::template Subspace<1>;

        PhaseFieldForBrittleFractures(FunctionSpace &space)
        : space_(space)
        {}

        void assemble(
            Vector &x,
            Matrix &H,
            Vector &g,
            Scalar &val
        )
        {
            USpace U;
            space_.subspace(0, U);
            CSpace C = space_.subspace(Dim);

            if(empty(H)) {
                space_.create_matrix(H);
            } else {
                H *= 0.0;
            }

            if(empty(g)) {
                space_.create_vector(g);
            } else {
                g.set(0.0);
            }

            FEFunction<FunctionSpace> x_fun(space_, x);
            auto x_coeff = x_fun.coefficient();

            val = 0.0;

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto x_view = x_coeff.view_device();

                Device::parallel_reduce(
                    space_.local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &)
                    {
                        return 0.0;
                    },
                    val
                );
            }

            val = x.comm().sum(val);
        }

    private:
        FunctionSpace space_;
    };

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
        using MatrixDxD      = utopia::StaticMatrix<Scalar, Dim, Dim>;

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
                MatrixDxD strain_n;
                MatrixDxD strain_p;

                ElemView e;
                space_view.elem(i, e);

                auto el_strain = strain_view.make(e);
                auto dx        = differential_view.make(e);

                const SizeType n_qp = el_strain.values.size();

                e_energy.set(0.0);

                for(SizeType qp = 0; qp < n_qp; ++qp) {
                    //reset strain tensor
                    strain_n.set(0.0);
                    strain_p.set(0.0);

                    //compute splitted strain
                    for(int d = 0; d < Dim; ++d) {
                        auto e_val = el_strain.values[qp][d];
                        el_strain.vectors[qp].col(d, v);

                        auto outer_v = outer(v, v);

                        auto eig_p = split_p(e_val);
                        auto eig_n = split_m(e_val);

                        strain_n += eig_n * outer_v;
                        strain_p += eig_p * outer_v;
                    }

                    const Scalar trace_e_n = trace(strain_n);
                    const Scalar trace_e_p = trace(strain_p);

                    e_energy[0] += ((0.5 * lambda * trace_e_n * trace_e_n) +
                               (mu * inner(strain_n, strain_n))) * dx(qp);

                    e_energy[1] += ((0.5 * lambda * trace_e_p * trace_e_p) +
                               (mu * inner(strain_p, strain_p))) * dx(qp);
                }

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