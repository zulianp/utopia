#ifndef UTOPIA_INTREPID2_GRADIENT_HPP
#define UTOPIA_INTREPID2_GRADIENT_HPP

#include "utopia_Traits.hpp"
#include "utopia_Views.hpp"

#include "utopia_intrepid2_Commons.hpp"
#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_Field.hpp"

namespace utopia {
    namespace intrepid2 {

        template <typename Scalar>
        class Gradient : public Field<Scalar> {
        public:
            using Super = utopia::intrepid2::Field<Scalar>;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using ExecutionSpace = typename FE::ExecutionSpace;

            Gradient(const std::shared_ptr<FE> &fe, const std::string &name = "Gradient") : Super(fe) {
                this->set_name(name);
                this->set_tensor_size(fe->spatial_dimension());
            }

            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Op(const DynRankView &grad, const DynRankView &coeff)
                    : grad(grad), coeff(coeff), num_fields(grad.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp, const int d) const {
                    Scalar ret = 0.0;
                    for (int i = 0; i < num_fields; ++i) {
                        ret += coeff(cell, i) * grad(cell, i, qp, d);
                    }

                    return ret;
                }

                const DynRankView grad, coeff;
                const int num_fields;
            };

            class OpAndStore {
            public:
                UTOPIA_INLINE_FUNCTION OpAndStore(const DynRankView &grad, const DynRankView &coeff, DynRankView &field)
                    : op_(grad, coeff), field(field) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int qp, const int d) const {
                    field(cell, qp, d) = op_(cell, qp, d);
                }

                Op op_;
                DynRankView field;
            };

            void init(const Field<Scalar> &coeff) {
                assert(coeff.tensor_size() == 1);
                ensure_field();
                init(coeff.data());
            }

            void init(const DynRankView &coeff) {
                ensure_field();
                Kokkos::parallel_for(
                    this->name() + "::init", grad_range(), OpAndStore(this->fe()->grad, coeff, this->data()));
            }

            void ensure_field() override {
                if (this->data().extent(0) < this->fe()->num_cells() || this->data().extent(1) < this->fe()->num_qp() ||
                    this->data().extent(2) < this->fe()->spatial_dimension()) {
                    this->data() = DynRankView(
                        this->name(), this->fe()->num_cells(), this->fe()->num_qp(), this->fe()->spatial_dimension());
                } else {
                    fill(this->data(), 0.0);
                }
            }

            inline Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace> grad_range() {
                int num_cells = this->fe()->num_cells();
                int num_qp = this->fe()->num_qp();
                int spatial_dimension = this->fe()->spatial_dimension();

                return Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0},
                                                                              {num_cells, num_qp, spatial_dimension});
            }
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_GRADIENT_HPP
