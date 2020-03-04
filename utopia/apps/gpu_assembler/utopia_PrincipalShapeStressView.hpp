#ifndef UTOPIA_PRINCIPAL_SHAPE_STRESS_VIEW_HPP
#define UTOPIA_PRINCIPAL_SHAPE_STRESS_VIEW_HPP

#include "utopia_LaplacianView.hpp"
#include "utopia_Utils.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_split_matrix.hpp"
#include "utopia_DeviceIdentity.hpp"

namespace utopia {


    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PrincipalShapeStress {};

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class PrincipalShapeStress< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:

        using FunctionSpace           = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector                  = typename FunctionSpace::Vector;
        using Scalar                  = typename FunctionSpace::Scalar;

        using FunctionSpaceViewDevice   = typename FunctionSpace::ViewDevice;


        using Elem      = typename FunctionSpace::ViewDevice::Elem;
        using GradValue = typename Elem::GradValue;
        static const int Dim = Elem::Dim;

        static const int NFunctions = Elem::NFunctions;

        class ViewDevice {
        public:
            // ArrayView<GradValue, NFunctions, Quadrature::NPoints> positive, negative;
            ArrayView<GradValue, NFunctions, Quadrature::NPoints> stress;
            ViewDevice(){}

            ViewDevice(const ViewDevice &other)
            {
                for(int j = 0; j < NFunctions; ++j) {
                    for(int i = 0; i < Quadrature::NPoints; ++i) {
                        stress(j, i).copy(other.stress(j, i));
                        // positive(j, i).copy(other.positive(j, i));
                    }
                }
            }
        };

        PrincipalShapeStress(
            const FunctionSpace &space,
            const Quadrature &q,
            const Scalar &mu,
            const Scalar &lambda)

        {
            init(space, q, mu, lambda);
        }

        inline const ViewDevice &view_device() const
        {
            return view_device_;
        }

    private:
        ViewDevice view_device_;

        void compute_aggregate_stress(
            const FunctionSpace &space,
            const Quadrature &q,
            const Scalar &mu,
            const Scalar &lambda //,
            // ArrayView<GradValue, NFunctions, Quadrature::NPoints> &stress
        )
        {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);
            auto grad_view = grad.view_host();

            Elem e;
            space.elem(0, e);

            auto g  = grad_view.make(e);

            GradValue strain;

            for(SizeType i = 0; i < e.n_functions(); ++i) {
                for(SizeType k = 0; k < Quadrature::NPoints; ++k) {
                    g.get(i, k, strain);
                    strain.symmetrize();

                    view_device_.stress(i, k) = 2.0 * mu * strain + lambda * trace(strain) * (device::identity<Scalar>());
                }
            }
        }

        // bool check_stress(const FunctionSpace &space, const Quadrature &q, const Scalar &mu, const Scalar &lambda)
        // {
        //     ArrayView<GradValue, NFunctions, Quadrature::NPoints> stress_qp;
        //     compute_aggregate_stress(space, q, mu, lambda, stress_qp);

        //     bool ok = true;
        //     GradValue actual_stress;
        //     for(SizeType i = 0; i < NFunctions; ++i) {
        //         for(SizeType k = 0; k < Quadrature::NPoints; ++k) {
        //             actual_stress = view_device_.positive(i, k) - view_device_.negative(i, k);
        //             if(!approxeq(stress_qp(i, k), actual_stress, 1e-8)) {
        //                 disp(std::to_string(i) + "," + std::to_string(k) + ")");
        //                 disp(stress_qp(i, k));
        //                 disp("==");
        //                 disp(actual_stress);
        //                 disp("neg");
        //                 disp(view_device_.negative(i, k));
        //                 disp("pos");
        //                 disp(view_device_.positive(i, k));
        //                 ok = false;
        //                 UTOPIA_DEVICE_ASSERT(false);
        //             }
        //         }
        //     }

        //     UTOPIA_DEVICE_ASSERT(ok);
        //     return ok;
        // }

        void init(const FunctionSpace &space, const Quadrature &q, const Scalar &mu, const Scalar &lambda)
        {
            compute_aggregate_stress(space, q, mu, lambda);
        }

        // void init(const FunctionSpace &space, const Quadrature &q, const Scalar &mu, const Scalar &lambda)
        // {
        //     PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);
        //     auto grad_view = grad.view_host();

        //     Elem e;
        //     space.elem(0, e);

        //     auto g  = grad_view.make(e);

        //     GradValue strain, strain_p, strain_n;

        //     StaticVector<Scalar, Dim> values;
        //     StaticMatrix<Scalar, Dim, Dim> vectors;
        //     StaticVector<Scalar, Dim> v;

        //     for(SizeType k = 0; k < Quadrature::NPoints; ++k) {
        //         for(SizeType i = 0; i < e.n_functions(); ++i) {
        //             g.get(i, k, strain);
        //             strain.symmetrize();

        //             // if(i == 0 && k == 2) {
        //             //     std::cout << "EHRERE" << std::endl;
        //             // }

        //             eig(strain, values, vectors);



        //             const Scalar sum_eigs = sum(values);
        //             const Scalar sum_eigs_p = ramp_fun_positive(sum_eigs);
        //             const Scalar sum_eigs_n = ramp_fun_negative(sum_eigs);

        //             // if(i == 0 && k == 2) {
        //             //     disp("-------");
        //             //     disp(strain);
        //             //     disp("eig: ");
        //             //     disp(values);
        //             //     disp(vectors);
        //             //     disp("-------");

        //             //     disp(sum_eigs, "sum_eigs");
        //             //     disp(sum_eigs_n, "sum-");
        //             //     disp(sum_eigs_p, "sum+");
        //             //     disp("------");
        //             // }

        //             auto &stress_p = view_device_.positive(i, k);
        //             auto &stress_n = view_device_.negative(i, k);

        //             stress_p.set(0.0);
        //             stress_n.set(0.0);

        //             for(int d = 0; d < Dim; ++d) {
        //                 vectors.col(d, v);
        //                 auto outer_v = outer(v, v);

        //                 const Scalar eig_p = ramp_fun_positive(values[d]);
        //                 const Scalar eig_n = ramp_fun_negative(values[d]);

        //                 // if(k == 0 && i == 0) {
        //                 //     disp(eig_n);
        //                 // }

        //                 stress_n += (lambda * sum_eigs_n + 2.0 * mu * eig_n) * outer_v;
        //                 stress_p += (lambda * sum_eigs_p + 2.0 * mu * eig_p) * outer_v;

        //             }
        //         }
        //     }

        //     UTOPIA_DEVICE_ASSERT( check_stress(space, q, mu, lambda) );
        // }
    };

}

#endif //UTOPIA_PRINCIPAL_SHAPE_STRESS_VIEW_HPP