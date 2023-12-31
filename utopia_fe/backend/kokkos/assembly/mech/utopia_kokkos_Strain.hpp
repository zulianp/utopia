#ifndef UTOPIA_KOKKOS_STRAIN_HPP
#define UTOPIA_KOKKOS_STRAIN_HPP

#include "utopia_Traits.hpp"

#include "utopia_kokkos_Gradient.hpp"
#include "utopia_kokkos_StrainOp.hpp"

namespace utopia {
    namespace kokkos {

        template <typename FE_>
        class Strain : public QPTensorField<FE_> {
        public:
            using Super = utopia::kokkos::QPTensorField<FE_>;
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using ExecutionSpace = typename FE::ExecutionSpace;

            bool is_coefficient() const override { return false; }

            void init_linearized(const Field<FE> &coeff) {
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
                // using Op = typename LinearizedStrain<Scalar, Dim>::Interpolate;
                using Op = utopia::kokkos::kernels::
                    InterpolateLinearizedStrainOp<Dim, Scalar, typename FE::Gradient, DynRankView>;

                using StoreOp = utopia::kokkos::kernels::StoreInterpolatedStrain<Op, DynRankView>;

                UTOPIA_INLINE_FUNCTION InterpolateAndStoreLinearized(const DynRankView &grad,
                                                                     const DynRankView &coeff,
                                                                     DynRankView &data)
                    : op(Op(grad, coeff), data) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                       const int qp,
                                                       const int sub_i,
                                                       const int sub_j) const {
                    op(cell, qp, sub_i, sub_j);
                }

                StoreOp op;
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

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_STRAIN_HPP
