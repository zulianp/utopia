#ifndef UTOPIA_DEVICE_EIGEN_DECOMPOSITION_NEW_HPP
#define UTOPIA_DEVICE_EIGEN_DECOMPOSITION_NEW_HPP

#include "utopia_DeviceFind.hpp"
#include "utopia_DeviceIdentity.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

namespace utopia {

    template <class Expr>
    class DeviceEigenDecompositionNew : public DeviceExpression<DeviceEigenDecompositionNew<Expr>> {
    public:
        using Scalar = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;
        using Matrix2x2 = utopia::StaticMatrix2x2<Scalar>;
        using Matrix3x3 = utopia::StaticMatrix3x3<Scalar>;
        using Vector3 = utopia::StaticVector3<Scalar>;

        template <class Vector, class ResultMat>
        UTOPIA_INLINE_FUNCTION static void apply(const Expr &mat, Vector &eigen_values, ResultMat &eigen_vectors) {
            switch (eigen_values.size()) {
                case 2: {
                    apply_2(mat, eigen_values, eigen_vectors);
                    return;
                }
                case 3: {
                    apply_3(mat, eigen_values, eigen_vectors);
                    return;
                }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return;
                }
            }

            UTOPIA_DEVICE_ASSERT(check(mat, eigen_values, eigen_vectors));
        }

    private:
        template <class MatT, class VecT>
        static UTOPIA_INLINE_FUNCTION void eigenvalues_2(const MatT &m, VecT &roots) {
            const Scalar t0 =
                Scalar(0.5) * device::sqrt(device::squared(m(0, 0) - m(1, 1)) + Scalar(4) * device::squared(m(1, 0)));
            const Scalar t1 = Scalar(0.5) * (m(0, 0) + m(1, 1));
            roots(0) = t1 - t0;
            roots(1) = t1 + t0;
        }

        template <class EigenValues, class EigenVectors>
        UTOPIA_INLINE_FUNCTION static void apply_2(const Expr &mat,
                                                   EigenValues &eigen_values,
                                                   EigenVectors &eigen_vectors) {
            Scalar scale = max(abs(mat));
            scale = device::max(scale, Scalar(1));
            Matrix2x2 scaled_mat = mat / scale;

            // Compute the eigenvalues
            eigenvalues_2(scaled_mat, eigen_values);

            if ((eigen_values(1) - eigen_values(0)) <= device::abs(eigen_values(1)) * device::epsilon<Scalar>()) {
                eigen_vectors.identity();
            } else {
                scaled_mat -= eigen_values(1) * device::identity<Scalar>();
                Scalar a2 = device::squared(scaled_mat(0, 0));
                Scalar c2 = device::squared(scaled_mat(1, 1));
                Scalar b2 = device::squared(scaled_mat(1, 0));

                if (a2 > c2) {
                    const Scalar norm_col = device::sqrt(a2 + b2);
                    UTOPIA_DEVICE_ASSERT(norm_col > 0);
                    eigen_vectors(0, 1) = -scaled_mat(1, 0) / norm_col;
                    eigen_vectors(1, 1) = scaled_mat(0, 0) / norm_col;
                } else {
                    const Scalar norm_col = device::sqrt(c2 + b2);
                    UTOPIA_DEVICE_ASSERT(norm_col > 0);
                    eigen_vectors(0, 1) = -scaled_mat(1, 1) / norm_col;
                    eigen_vectors(1, 1) = scaled_mat(1, 0) / norm_col;
                }

                //((-y, x))
                eigen_vectors(0, 0) = -eigen_vectors(1, 1);
                eigen_vectors(1, 0) = eigen_vectors(0, 1);
            }

            eigen_values *= scale;
        }

        //////////////////////////////////////////////////////////////

        template <class MatT, class VecT>
        static UTOPIA_INLINE_FUNCTION void eigenvalues_3(const MatT &m, VecT &roots) {
            const Scalar s_inv3 = Scalar(1.0) / Scalar(3.0);
            const Scalar s_sqrt3 = device::sqrt(Scalar(3.0));

            // The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
            // eigenvalues are the roots to this equation, all guaranteed to be
            // real-valued, because the matrix is symmetric.
            Scalar c0 = m(0, 0) * m(1, 1) * m(2, 2) + Scalar(2) * m(1, 0) * m(2, 0) * m(2, 1) -
                        m(0, 0) * m(2, 1) * m(2, 1) - m(1, 1) * m(2, 0) * m(2, 0) - m(2, 2) * m(1, 0) * m(1, 0);
            Scalar c1 = m(0, 0) * m(1, 1) - m(1, 0) * m(1, 0) + m(0, 0) * m(2, 2) - m(2, 0) * m(2, 0) +
                        m(1, 1) * m(2, 2) - m(2, 1) * m(2, 1);
            Scalar c2 = m(0, 0) + m(1, 1) + m(2, 2);

            // Construct the parameters used in classifying the roots of the equation
            // and in solving the equation for the roots in closed form.
            Scalar c2_over_3 = c2 * s_inv3;
            Scalar a_over_3 = (c2 * c2_over_3 - c1) * s_inv3;
            if (a_over_3 < Scalar(0)) a_over_3 = Scalar(0);

            Scalar half_b = Scalar(0.5) * (c0 + c2_over_3 * (Scalar(2) * c2_over_3 * c2_over_3 - c1));

            Scalar q = a_over_3 * a_over_3 * a_over_3 - half_b * half_b;
            if (q < Scalar(0)) q = Scalar(0);

            // Compute the eigenvalues by solving for the roots of the polynomial.
            Scalar rho = device::sqrt(a_over_3);
            Scalar theta = device::atan2(device::sqrt(q), half_b) *
                           s_inv3;  // since sqrt(q) > 0, atan2 is in [0, pi] and theta is in [0, pi/3]
            Scalar cos_theta = device::cos(theta);
            Scalar sin_theta = device::sin(theta);
            // roots are already sorted, since cos is monotonically decreasing on [0, pi]
            roots(0) = c2_over_3 - rho * (cos_theta + s_sqrt3 * sin_theta);  // == 2*rho*cos(theta+2pi/3)
            roots(1) = c2_over_3 - rho * (cos_theta - s_sqrt3 * sin_theta);  // == 2*rho*cos(theta+ pi/3)
            roots(2) = c2_over_3 + Scalar(2) * rho * cos_theta;
        }

        template <class MatT, class VecT>
        static inline void extract_kernel_3(const MatT &mat, VecT &res, VecT &representative) {
            SizeType i0 = imax(abs(diag(mat)));
            mat.col(i0, representative);

            Scalar n0, n1;
            Vector3 temp, c0, c1;

            mat.col((i0 + 1) % 3, temp);
            c0 = cross(representative, temp);

            mat.col((i0 + 2) % 3, temp);
            c1 = cross(representative, temp);

            n0 = dot(c0, c0);
            n1 = dot(c1, c1);

            if (n0 > n1)
                res = c0 / device::sqrt(n0);
            else
                res = c1 / device::sqrt(n1);
        }

        template <class EigenValues, class EigenVectors>
        UTOPIA_INLINE_FUNCTION static void apply_3(const Expr &mat,
                                                   EigenValues &eigen_values,
                                                   EigenVectors &eigen_vectors) {
            Vector3 vk, vl;
            Matrix3x3 scaled_mat;
            Matrix3x3 tmp;

            // Shift the matrix to the mean eigenvalue and map the matrix coefficients to [-1:1] to avoid over- and
            // underflow.
            Scalar shift = trace(mat) / Scalar(3);

            scaled_mat.copy(mat);
            scaled_mat -= shift * device::identity<Scalar>();

            Scalar scale = max(abs(scaled_mat));
            if (scale > 0) scaled_mat /= scale;

            eigenvalues_3(scaled_mat, eigen_values);

            if ((eigen_values(2) - eigen_values(0)) <= device::epsilon<Scalar>()) {
                eigen_vectors.identity();
            } else {
                tmp.copy(scaled_mat);

                // Compute the eigenvector of the most distinct eigenvalue
                Scalar d0 = eigen_values(2) - eigen_values(1);
                Scalar d1 = eigen_values(1) - eigen_values(0);
                SizeType k(0), l(2);
                if (d0 > d1) {
                    device::swap(k, l);
                    d0 = d1;
                }

                // Compute the eigenvector of index k
                {
                    tmp -= eigen_values(k) * device::identity<Scalar>();
                    extract_kernel_3(tmp, vk, vl);
                    eigen_vectors.set_col(k, vk);
                    eigen_vectors.set_col(l, vl);
                }

                // Compute eigenvector of index l
                if (d0 <= 2 * device::epsilon<Scalar>() * d1) {
                    // If d0 is too small, then the two other eigenvalues are numerically the same,
                    // and thus we only have to ortho-normalize the near orthogonal vector we saved above.
                    vl -= dot(vk, vl) * vl;
                    vl /= norm2(vl);
                    eigen_vectors.set_col(l, vl);
                } else {
                    tmp.copy(scaled_mat);
                    tmp -= eigen_values(l) * device::identity<Scalar>();
                    extract_kernel_3(tmp, vl, vk);
                    eigen_vectors.set_col(l, vl);
                }

                eigen_vectors.col(0, vk);
                eigen_vectors.col(2, vl);

                Vector3 v1 = cross(vl, vk);
                v1 /= norm2(v1);
                eigen_vectors.set_col(1, v1);
            }

            // Return to orginal frame
            eigen_values *= scale;
            eigen_values += shift;
        }

        template <class Mat, class Values, class Vectors>
        UTOPIA_INLINE_FUNCTION static bool check(const Mat &mat, const Values &values, const Vectors &vectors) {
            auto sum_v = sum(vectors);
            auto sum_e = sum(values);

            bool ok = !device::isnan(sum_v) && !device::isnan(sum_e);

            UTOPIA_DEVICE_ASSERT(!device::isnan(sum_v));
            UTOPIA_DEVICE_ASSERT(!device::isnan(sum_e));

            if (ok) {
                // check recomposition
                ok = approxeq(mat, vectors * diag(values) * transpose(vectors), 1e-8);
            }

            // if (!ok) {
            //     disp("--------------------");
            //     disp(mat);
            //     disp("--------------------");
            //     disp(vectors);
            //     disp("--------------------");
            //     disp(values);
            //     disp("--------------------");

            //     Mat reconstructed = vectors * diag(values) * transpose(vectors);
            //     disp(reconstructed);
            // }

            UTOPIA_DEVICE_ASSERT(ok);
            return ok;
        }
    };

}  // namespace utopia
#endif  // UTOPIA_DEVICE_EIGEN_DECOMPOSITION_NEW_HPP
