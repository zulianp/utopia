#ifndef UTOPIA_MARS_UNIFORM_FE_HPP
#define UTOPIA_MARS_UNIFORM_FE_HPP

#include "mars.hpp"
#include "utopia_kokkos_Traits.hpp"
#include "utopia_mars_Hex8.hpp"
#include "utopia_mars_Hex8Quadrature.hpp"

#include "utopia_kokkos_UniformFE.hpp"

namespace utopia {

    template <int Dim, typename Scalar_>
    class Traits<Scalar_[Dim]> {
    public:
        using Scalar = Scalar_;
    };

    namespace mars {

        template <typename Scalar, ::mars::Integer XDim_, ::mars::Integer YDim_, ::mars::Integer ZDim_>
        using Texture3D = Kokkos::View<Scalar[XDim_][YDim_][ZDim_],
                                       ::mars::KokkosLayout,
                                       ::mars::KokkosSpace,
                                       Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

        template <class ElemType, class Quadrature>
        class UniformFE {};

        template <typename Scalar, int PhysicalDim, class Quadrature>
        class UniformFE<::mars::FEQuad4<Scalar, PhysicalDim>, Quadrature> {
        public:
            using Elem = ::mars::FEQuad4<Scalar, PhysicalDim>;

            static constexpr int ManifoldDim = 2;
            static constexpr int NQPoints = Elem::n_qp;
            static constexpr int NFun = 4;

            using Coordinates = ::mars::ViewMatrixTextureC<Scalar, 4, PhysicalDim>;

            ::mars::ViewMatrixTextureC<Scalar, ManifoldDim, PhysicalDim> J_inv;
            ::mars::ViewVectorTextureC<Scalar, NQPoints> measure;

            Coordinates coordinates;
            ::mars::ViewObject<Scalar> det_J;

            ::mars::ViewMatrixTextureC<Scalar, 4, NQPoints> fun;
            Texture3D<Scalar, 4, NQPoints, PhysicalDim> grad;

            UniformFE() : J_inv("J_inv"), measure("measure"), det_J("det_J"), fun("fun"), grad("grad") {}

            void init(const Coordinates &in_coordinates) {
                coordinates = in_coordinates;
                auto determinant_j = det_J;
                auto inverse_j = J_inv;
                auto meas = measure;

                Kokkos::parallel_for(
                    1, MARS_LAMBDA(const int) {
                        Scalar J_e[PhysicalDim * ManifoldDim], J_inv_e[ManifoldDim * PhysicalDim];

                        // col 0, p1
                        for (int d = 0; d < PhysicalDim; ++d) {
                            J_e[d * ManifoldDim] = in_coordinates(1, d) - in_coordinates(0, d);
                        }

                        // col 1, p3
                        for (int d = 0; d < PhysicalDim; ++d) {
                            J_e[d * ManifoldDim + 1] = in_coordinates(3, d) - in_coordinates(0, d);
                        }

                        Scalar det_J_e;
                        ::mars::Invert<ManifoldDim>::apply(J_e, J_inv_e, det_J_e);

                        for (int d1 = 0; d1 < ManifoldDim; ++d1) {
                            for (int d2 = 0; d2 < PhysicalDim; ++d2) {
                                inverse_j(d1, d2) = J_inv_e[d1 * PhysicalDim + d2];
                            }
                        }

                        assert(det_J_e > 0);
                        determinant_j(0) = det_J_e;
                    });

                Quadrature q = Quadrature::make();
                auto points = q.q_p;
                auto weights = q.q_w;

                auto fun = this->fun;
                auto grad = this->grad;

                Kokkos::parallel_for(
                    NQPoints, MARS_LAMBDA(const int qp) {
                        Scalar J_inv_e[ManifoldDim * PhysicalDim];
                        Scalar p_qp[ManifoldDim], grad_physical[ManifoldDim];

                        assert(determinant_j(0) > 0);
                        assert(weights(qp) != 0.0);

                        meas(qp) = determinant_j(0) * weights(qp);

                        for (int d1 = 0; d1 < ManifoldDim; ++d1) {
                            p_qp[d1] = points(qp, d1);

                            for (int d2 = 0; d2 < PhysicalDim; ++d2) {
                                J_inv_e[d1 * PhysicalDim + d2] = inverse_j(d1, d2);
                            }
                        }

                        for (int i = 0; i < NFun; ++i) {
                            fun(i, qp) = Elem::Fun::f(i, p_qp);

                            Elem::Grad::affine_f(i, J_inv_e, p_qp, grad_physical);

                            for (int d = 0; d < PhysicalDim; ++d) {
                                grad(i, qp, d) = grad_physical[d];
                            }
                        }
                    });
            }
        };

        template <typename Scalar, class Quadrature>
        class UniformFE<UniformHex8<Scalar>, Quadrature> {
        public:
            using Elem = utopia::UniformHex8<Scalar>;

            static constexpr int PhysicalDim = 3;
            static constexpr int ManifoldDim = 3;
            static constexpr int NQPoints = Quadrature::NPoints;
            static constexpr int NFun = 8;

            using Coordinates = ::mars::ViewMatrixTextureC<Scalar, 8, PhysicalDim>;
            ::mars::ViewMatrixTextureC<Scalar, ManifoldDim, PhysicalDim> J_inv;
            ::mars::ViewVectorTextureC<Scalar, NQPoints> measure;

            Coordinates coordinates;
            ::mars::ViewObject<Scalar> det_J;

            ::mars::ViewMatrixTextureC<Scalar, 8, NQPoints> fun;
            Texture3D<Scalar, 8, NQPoints, PhysicalDim> grad;

            UniformFE() : J_inv("J_inv"), measure("measure"), det_J("det_J"), fun("fun"), grad("grad") {}

            void init(const Coordinates &in_coordinates) {
                coordinates = in_coordinates;
                auto determinant_j = det_J;
                auto inverse_j = J_inv;
                auto meas = measure;

                ::mars::ViewMatrixTextureC<Scalar, NQPoints, PhysicalDim> points("q_points");
                ::mars::ViewVectorTextureC<Scalar, NQPoints> weights("q_weights");

                Kokkos::parallel_for(
                    1, MARS_LAMBDA(const int) {
                        Scalar J_e[PhysicalDim * ManifoldDim], J_inv_e[ManifoldDim * PhysicalDim];

                        // col 0, p1
                        for (int d = 0; d < PhysicalDim; ++d) {
                            J_e[d * ManifoldDim] = in_coordinates(1, d) - in_coordinates(0, d);
                        }

                        // col 1, p3
                        for (int d = 0; d < PhysicalDim; ++d) {
                            J_e[d * ManifoldDim + 1] = in_coordinates(3, d) - in_coordinates(0, d);
                        }

                        // col 2, p4
                        for (int d = 0; d < PhysicalDim; ++d) {
                            J_e[d * ManifoldDim + 2] = in_coordinates(4, d) - in_coordinates(0, d);
                        }

                        Scalar det_J_e;
                        ::mars::Invert<ManifoldDim>::apply(J_e, J_inv_e, det_J_e);

                        for (int d1 = 0; d1 < ManifoldDim; ++d1) {
                            for (int d2 = 0; d2 < PhysicalDim; ++d2) {
                                inverse_j(d1, d2) = J_inv_e[d1 * PhysicalDim + d2];
                                assert(inverse_j(d1, d2) == inverse_j(d1, d2));
                            }

                            assert(inverse_j(d1, d1) != 0);
                        }

                        assert(det_J_e == det_J_e);
                        assert(det_J_e > 0);
                        determinant_j(0) = det_J_e;

                        Quadrature::get(points, weights);
                    });

                auto fun = this->fun;
                auto grad = this->grad;

                Kokkos::parallel_for(
                    NQPoints, MARS_LAMBDA(const int qp) {
                        meas(qp) = determinant_j(0) * weights(qp);

                        Elem elem;
                        StaticVector<Scalar, PhysicalDim> h, translation;
                        StaticVector<Scalar, ManifoldDim> p_qp;
                        StaticVector<Scalar, PhysicalDim> grad_physical;

                        for (int d2 = 0; d2 < PhysicalDim; ++d2) {
                            assert(inverse_j(d2, d2) != 0);
                            h[d2] = 1. / inverse_j(d2, d2);
                            translation[d2] = in_coordinates(0, d2);
                            assert(h[d2] == h[d2]);
                        }

                        elem.set(translation, h);

                        for (int d1 = 0; d1 < ManifoldDim; ++d1) {
                            p_qp[d1] = points(qp, d1);
                        }

                        for (int i = 0; i < NFun; ++i) {
                            fun(i, qp) = Elem::fun(i, p_qp);

                            elem.grad(i, p_qp, grad_physical);

                            for (int d = 0; d < PhysicalDim; ++d) {
                                grad(i, qp, d) = grad_physical[d];
                                assert(grad(i, qp, d) == grad(i, qp, d));
                                // assert(std::isfinite(grad_physical[d]));
                            }
                        }
                    });
            }
        };

        template <typename Scalar, int Dim>
        class FETypeSelect {};

        template <typename Scalar>
        class FETypeSelect<Scalar, 2> {
        public:
            using Elem = ::mars::FEQuad4<Scalar, 2>;
            using Quadrature = typename Elem::Quadrature;
            using Type = utopia::mars::UniformFE<Elem, Quadrature>;
        };

        template <typename Scalar>
        class FETypeSelect<Scalar, 3> {
        public:
            using Elem = utopia::UniformHex8<Scalar>;
            using Quadrature = utopia::Hex8Quadrature<Scalar, 2, 3, 27>;
            using Type = utopia::mars::UniformFE<Elem, Quadrature>;
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_UNIFORM_FE_HPP
