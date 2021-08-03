#ifndef UTOPIA_KOKKOS_TRANSPORT_OP_HPP
#define UTOPIA_KOKKOS_TRANSPORT_OP_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <int Dim, typename Scalar, class Field, class Grad, class Fun, class Measure>
            class TransportOp {
            public:
                UTOPIA_INLINE_FUNCTION TransportOp(const Field &vector_field,
                                                   const Grad &grad,
                                                   const Fun &fun,
                                                   const Measure &measure)
                    : vector_field(vector_field), grad(grad), fun(fun), measure(measure), n_qp(measure.extent(1)) {
                    assert(vector_field.extent(0) == measure.extent(0));
                    assert(vector_field.extent(1) == measure.extent(1));
                    assert(vector_field.extent(2) == grad.extent(3));
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i, const int &j) const {
                    Scalar integral = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar val = 0.0;
                        for (int dj = 0; dj < Dim; ++dj) {
                            assert(vector_field(cell, qp, dj) == vector_field(cell, qp, dj));
                            assert(grad(cell, j, qp, dj) == grad(cell, j, qp, dj));

                            val += grad(cell, j, qp, dj) * vector_field(cell, qp, dj);
                        }

                        assert(val == val);

                        auto dX = measure(cell, qp);
                        assert(dX == dX);

                        val *= fun(i, qp) * dX;
                        integral += val;
                    }

                    assert(integral == integral);
                    return integral;
                }

                Field vector_field;
                Grad grad;
                Fun fun;
                Measure measure;
                const int n_qp;
            };

        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_TRANSPORT_OP_HPP
