#ifndef UTOPIA_KOKKOS_PRINCIPAL_STRESSES_HPP
#define UTOPIA_KOKKOS_PRINCIPAL_STRESSES_HPP

namespace utopia {
    namespace kokkos {

        template <int Dim, typename Scalar, class StressKernel, class Measure, class Output>
        class StorePrincipalStress {
        public:
            StorePrincipalStress(const StressKernel &stress, const Measure &measure, Output &output)
                : stress(stress), measure(measure), output(output), n_quad_points(measure.extent(1)) {}

            UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int qp) const {
                StaticMatrix<Scalar, Dim, Dim> mat;
                StaticVector<Scalar, Dim> e;

                for (int d1 = 0; d1 < Dim; ++d1) {
                    for (int d2 = 0; d2 < Dim; ++d2) {
                        mat(d1, d2) = stress(cell, qp, d1, d2);
                    }
                }

                utopia::eig(mat, e);

                auto dX = measure(cell, qp);
                for (int d = 0; d < Dim; ++d) {
                    output(cell, d) += e[d] * dX;
                }
            }

            UTOPIA_INLINE_FUNCTION void operator()(const int cell) const {
                Scalar cell_measure = 0.0;
                for (int qp = 0; qp < n_quad_points; ++qp) {
                    (*this)(cell, qp);

                    cell_measure += measure(cell, qp);
                }

                for (int d = 0; d < Dim; ++d) {
                    output(cell, d) /= cell_measure;
                }
            }

            StressKernel stress;
            Measure measure;
            Output output;
            int n_quad_points;
        };
    }  // namespace kokkos

}  // namespace utopia

#endif  // UTOPIA_KOKKOS_PRINCIPAL_STRESSES_HPP