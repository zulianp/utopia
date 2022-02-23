#ifndef UTOPIA_KOKKOS_LINEAR_SIMPLEX_HPP
#define UTOPIA_KOKKOS_LINEAR_SIMPLEX_HPP

#include "utopia_Views.hpp"

namespace utopia {
    namespace kokkos {

        template <class View, int Dim>
        class Simplex {
        public:
            using SizeType = typename Traits<View>::SizeType;
            using Scalar = typename Traits<View>::Scalar;

            UTOPIA_INLINE_FUNCTION static Scalar fun0(const Scalar *local_p) {
                Scalar ret = 1;
                for (int d = 0; d < Dim; ++d) {
                    ret -= local_p[d];
                }

                return ret;
            }

            UTOPIA_INLINE_FUNCTION static Scalar fun(const int nf, const Scalar *local_p) {
                if (nf == 0) {
                    return fun0(local_p);
                }

                return local_p[nf - 1];
            }

            UTOPIA_INLINE_FUNCTION static Scalar fun(const int nf, const StaticVector<Scalar, Dim> &local_p) {
                if (nf == 0) {
                    return fun0(&local_p[0]);
                }

                return local_p[nf - 1];
            }

            UTOPIA_INLINE_FUNCTION void inverse_transform(const SizeType &cell,
                                                          const Scalar *global_p,
                                                          Scalar *local_p) const {
                StaticMatrix<Scalar, Dim, Dim> A, A_inv;
                StaticVector<Scalar, Dim> b, xg, xl;

                static constexpr int n_nodes = Dim + 1;
                assert(int(points.extent(1)) == n_nodes);

                for (int d = 0; d < Dim; ++d) {
                    xg[d] = global_p[d];
                }

                inverse_transform(cell, xg, xl);

                for (int d = 0; d < Dim; ++d) {
                    local_p[d] = xl[d];
                }
            }

            UTOPIA_INLINE_FUNCTION void inverse_transform(const SizeType &cell,
                                                          const StaticVector<Scalar, Dim> &xg,
                                                          StaticVector<Scalar, Dim> &xl) const {
                StaticMatrix<Scalar, Dim, Dim> A, A_inv;
                StaticVector<Scalar, Dim> b;

                static constexpr int n_nodes = Dim + 1;
                assert(int(points.extent(1)) == n_nodes);

                for (int d = 0; d < Dim; ++d) {
                    b[d] = points(cell, 0, d);
                }

                for (int k = 0; k < Dim; ++k) {
                    for (int d = 0; d < Dim; ++d) {
                        A(k, d) = points(cell, k + 1, d) - b[d];
                    }
                }

                A_inv = inv(A);

                b = xg - b;
                xl = A_inv * b;
            }

            UTOPIA_INLINE_FUNCTION Simplex(const View &points) : points(points) {}

            View points;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_LINEAR_SIMPLEX_HPP
