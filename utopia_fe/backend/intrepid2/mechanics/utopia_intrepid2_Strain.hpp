#ifndef UTOPIA_INTREPID2_STRAIN_HPP
#define UTOPIA_INTREPID2_STRAIN_HPP

#include "utopia_Traits.hpp"

#include "utopia_intrepid2_Gradient.hpp"

namespace utopia {
    namespace intrepid2 {

        template <typename Scalar, int Dim>
        class LinearizedStrain {
        public:
            template <class Grad>
            UTOPIA_INLINE_FUNCTION static constexpr auto inner(const Grad &grad,
                                                               const int cell,
                                                               const int i,
                                                               const int j,
                                                               const int qp,
                                                               const int sub_i,
                                                               const int sub_j) {
                if (sub_i == sub_j) {
                    Scalar ret = 0.0;
                    for (int d = 0; d < Dim; ++d) {
                        ret += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                    }

                    ret *= 0.5;
                    ret += 0.5 * grad(cell, i, qp, sub_i) * grad(cell, j, qp, sub_i);
                    return ret;

                } else {
                    Scalar ret = 0.5 * grad(cell, i, qp, sub_j) * grad(cell, j, qp, sub_i);
                    return ret;
                }
            }

            template <class Grad>
            UTOPIA_INLINE_FUNCTION static constexpr auto trace(const Grad &grad,
                                                               const int cell,
                                                               const int i,
                                                               const int qp,
                                                               const int sub_i) {
                return grad(cell, i, qp, sub_i);
            }

            class Interpolate {
            public:
                using GradOp = typename Gradient<Scalar>::Rank2Op;
                using FE = utopia::intrepid2::FE<Scalar>;
                using SizeType = typename FE::SizeType;
                using DynRankView = typename FE::DynRankView;

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
                                                         const int qp,
                                                         const int sub_i,
                                                         const int sub_j) const {
                    return 0.5 * (grad_(cell, qp, sub_i, sub_j) + grad_(cell, qp, sub_j, sub_i));
                }

                UTOPIA_INLINE_FUNCTION Scalar trace(const int cell, const int qp) const {
                    Scalar ret = (*this)(cell, qp, 0, 0);

                    for (int d = 1; d < Dim; ++d) {
                        ret += (*this)(cell, qp, d, d);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar squared_norm(const int cell, const int qp) const {
                    Scalar ret = 0.0;

                    for (int d1 = 0; d1 < Dim; ++d1) {
                        for (int d2 = 0; d2 < Dim; ++d2) {
                            auto x = (*this)(cell, qp, d1, d2);
                            ret += x * x;
                        }
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Interpolate(const DynRankView &grad, const DynRankView &coeff)
                    : grad_(grad, coeff) {}

                GradOp grad_;
            };

            // UTOPIA_INLINE_FUNCTION static void make(const int dim,
            //                                         Scalar *grad,
            //                                         StaticMatrix<Scalar, Dim, Dim> &strain) {
            //     strain.set(0.0);

            //     for (int i = 0; i < Dim; ++i) {
            //         strain(dim, i) = grad[i];
            //     }

            //     strain.symmetrize();
            // }
        };

        template <typename Scalar>
        class Strain : public QPTensorField<Scalar> {
        public:
            using Super = utopia::intrepid2::QPTensorField<Scalar>;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using ExecutionSpace = typename FE::ExecutionSpace;

            bool is_coefficient() const override { return false; }

            void init_linearized(const Field<Scalar> &coeff) {
                assert(coeff.tensor_size() > 0);
                assert(coeff.is_coefficient());
                this->set_tensor_size(coeff.tensor_size(), this->fe()->spatial_dimension());
                init_linearized(coeff.data());
            }

            Strain(const std::shared_ptr<FE> &fe, const std::string &name = "Strain") : Super(fe, name) {}
            // -----------
            // NVCC_PRIVATE:
            // compatibility with nvcc

            template <int Dim>
            class InterpolateAndStoreLinearized {
            public:
                using Op = typename LinearizedStrain<Scalar, 1>::Interpolate;

                UTOPIA_INLINE_FUNCTION InterpolateAndStoreLinearized(const DynRankView &grad,
                                                                     const DynRankView &coeff,
                                                                     DynRankView &data)
                    : op_(grad, coeff), data(data), n(grad.extent(3)) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                       const int qp,
                                                       const int sub_i,
                                                       const int sub_j) const {
                    data(cell, qp, sub_i * n + sub_j) = op_(cell, qp, sub_i, sub_j);
                }

                Op op_;
                DynRankView data;
                const int n;
            };

            void init_linearized(const DynRankView &coeff) {
                this->ensure_field();

                auto fe = this->fe();

                switch (fe->spatial_dimension()) {
                    case 1: {
                        Kokkos::parallel_for(this->name() + "::init_linearized",
                                             this->rank2_range(),
                                             InterpolateAndStoreLinearized<1>(this->fe()->grad(), coeff, this->data()));
                        break;
                    }

                    case 2: {
                        Kokkos::parallel_for(this->name() + "::init_linearized",
                                             this->rank2_range(),
                                             InterpolateAndStoreLinearized<2>(this->fe()->grad(), coeff, this->data()));
                        break;
                    }

                    case 3: {
                        Kokkos::parallel_for(this->name() + "::init_linearized",
                                             this->rank2_range(),
                                             InterpolateAndStoreLinearized<3>(this->fe()->grad(), coeff, this->data()));
                        break;
                    }

                    default: {
                        assert(false);
                        utopia::Utopia::Abort();
                    }
                }
            }
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_STRAIN_HPP
