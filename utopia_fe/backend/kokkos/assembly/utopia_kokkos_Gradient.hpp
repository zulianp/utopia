#ifndef UTOPIA_KOKKOS_GRADIENT_HPP
#define UTOPIA_KOKKOS_GRADIENT_HPP

#include "utopia_Traits.hpp"
#include "utopia_Views.hpp"

#include "utopia_kokkos_Commons.hpp"
#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_TensorField.hpp"

namespace utopia {
    namespace kokkos {

        template <typename FE_>
        class Gradient : public QPTensorField<FE_> {
        public:
            using FE = FE_;
            using Super = utopia::kokkos::QPTensorField<FE>;

            using Scalar = typename FE::Scalar;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using ExecutionSpace = typename FE::ExecutionSpace;

            Gradient(const std::shared_ptr<FE> &fe, const std::string &name = "Gradient") : Super(fe, name) {}

            class Rank1Op {
            public:
                UTOPIA_INLINE_FUNCTION Rank1Op(const DynRankView &grad, const DynRankView &coeff)
                    : grad(grad), coeff(coeff), n_shape_functions(grad.extent(1)), n_var(grad.extent(3)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp, const int d) const {
                    assert(d < n_var);

                    Scalar ret = 0.0;
                    for (int i = 0; i < n_shape_functions; ++i) {
                        ret += coeff(cell, i) * grad(cell, i, qp, d);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar squared_norm(const int cell, const int qp) const {
                    Scalar ret = 0.0;

                    for (int d = 0; d < n_var; ++d) {
                        auto x = (*this)(cell, qp, d);
                        ret += x * x;
                    }
                }

                const DynRankView grad, coeff;
                const int n_shape_functions;
                const int n_var;
            };

            class Rank1OpAndStore {
            public:
                UTOPIA_INLINE_FUNCTION Rank1OpAndStore(const DynRankView &grad,
                                                       const DynRankView &coeff,
                                                       DynRankView &field)
                    : op_(grad, coeff), field(field) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int qp, const int d) const {
                    field(cell, qp, d) = op_(cell, qp, d);
                }

                Rank1Op op_;
                DynRankView field;
            };

            class Rank2Op {
            public:
                UTOPIA_INLINE_FUNCTION Rank2Op(const DynRankView &grad, const DynRankView &coeff)
                    : grad(grad),
                      coeff(coeff),
                      n_shape_functions(grad.extent(1)),
                      n_var(coeff.extent(1) / n_shape_functions) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
                                                         const int qp,
                                                         const int var,
                                                         const int d) const {
                    Scalar ret = 0.0;
                    for (int i = 0; i < n_shape_functions; ++i) {
                        ret += coeff(cell, i * n_var + var) * grad(cell, i, qp, d);
                    }

                    return ret;
                }

                const DynRankView grad, coeff;
                const int n_shape_functions;
                const int n_var;
            };

            class Rank2OpAndStore {
            public:
                UTOPIA_INLINE_FUNCTION Rank2OpAndStore(const DynRankView &grad,
                                                       const DynRankView &coeff,
                                                       DynRankView &field)
                    : op_(grad, coeff), spatial_dim(grad.extent(3)), field(field) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int qp, const int var, const int d) const {
                    field(cell, qp, var * spatial_dim + d) = op_(cell, qp, var, d);
                }

                Rank2Op op_;
                const int spatial_dim;
                DynRankView field;
            };

            class Rank2SubOp {
            public:
                UTOPIA_INLINE_FUNCTION Rank2SubOp(const DynRankView &grad,
                                                  const DynRankView &coeff,
                                                  const int component)
                    : grad(grad),
                      coeff(coeff),
                      n_shape_functions(grad.extent(1)),
                      n_var(coeff.extent(1) / n_shape_functions),
                      component(component) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
                                                         const int qp,
                                                         const int var,
                                                         const int d) const {
                    Scalar ret = 0.0;
                    for (int i = 0; i < n_shape_functions; ++i) {
                        ret += coeff(cell, i * n_var + (component + var)) * grad(cell, i, qp, d);
                    }

                    return ret;
                }

                const DynRankView grad, coeff;
                const int n_shape_functions;
                const int n_var;
                const int component;
            };

            class Rank2SubOpAndStore {
            public:
                UTOPIA_INLINE_FUNCTION Rank2SubOpAndStore(const DynRankView &grad,
                                                          const DynRankView &coeff,
                                                          DynRankView &field,
                                                          const int component)
                    : op_(grad, coeff), spatial_dim(grad.extent(3)), field(field), component(component) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int qp, const int var, const int d) const {
                    field(cell, qp, var * spatial_dim + d) = op_(cell, qp, var, d);
                }

                Rank2Op op_;
                const int spatial_dim;
                DynRankView field;
                const int component;
            };

            void init(const Field<FE> &coeff) {
                assert(coeff.tensor_size() > 0);
                assert(coeff.is_coefficient());
                this->set_tensor_size(coeff.tensor_size(), this->fe()->spatial_dimension());
                init(coeff.data());
            }

            void init(const Field<FE> &coeff, int component, int size) {
                assert(coeff.tensor_size() > 0);
                assert(coeff.is_coefficient());

                this->set_tensor_size(size, this->fe()->spatial_dimension());
                init(coeff.data(), component);
            }

            bool is_coefficient() const override { return false; }

            // NVCC_PRIVATE:
            void init(const DynRankView &coeff) {
                this->ensure_field();

                if (this->rank() == 1) {
                    Kokkos::parallel_for(this->name() + "::init_rank1",
                                         this->rank1_range(),
                                         Rank1OpAndStore(this->fe()->grad(), coeff, this->data()));
                } else {
                    assert(this->rank() == 2);

                    Kokkos::parallel_for(this->name() + "::init_rank2",
                                         this->rank2_range(),
                                         Rank2OpAndStore(this->fe()->grad(), coeff, this->data()));
                }
            }

            void init(const DynRankView &coeff, int component) {
                this->ensure_field();

                if (this->rank() == 1) {
                    Kokkos::parallel_for(this->name() + "::init_rank1",
                                         this->rank1_range(),
                                         Rank1OpAndStore(this->fe()->grad(), coeff, this->data()));
                } else {
                    assert(this->rank() == 2);

                    Kokkos::parallel_for(this->name() + "::init_rank2",
                                         this->rank2_range(),
                                         Rank2SubOpAndStore(this->fe()->grad(), coeff, this->data(), component));
                }
            }
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_GRADIENT_HPP
