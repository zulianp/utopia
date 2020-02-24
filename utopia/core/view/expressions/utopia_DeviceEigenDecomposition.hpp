#ifndef UTOPIA_DEVICE_EIGEN_DECOMPOSITION_HPP
#define UTOPIA_DEVICE_EIGEN_DECOMPOSITION_HPP

#include "utopia_DeviceEigenValues.hpp"
#include "utopia_DeviceEigenVectors.hpp"

namespace utopia {
    //https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
    template<class Expr>
    class DeviceEigenDecomposition : public DeviceExpression<DeviceEigenDecomposition<Expr>> {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        template<class EigenValues, class EigenVectors>
        UTOPIA_INLINE_FUNCTION static void apply(
            const Expr &mat,
            EigenValues &eigen_values,
            EigenVectors &eigen_vectors)
        {

            Scalar tol = device::epsilon<Scalar>() * 1./max(abs(mat));
            if(mat.is_diagonal(tol)) {
                eigen_values = diag(mat);
                eigen_vectors.identity();
                return;
            }

            DeviceEigenValues<Expr>::apply(mat, eigen_values);
            DeviceEigenVectors<Expr>::apply(mat, eigen_values, eigen_vectors);

            UTOPIA_DEVICE_ASSERT(
                check(mat, eigen_values, eigen_vectors)
            );
        }

        template<class Mat, class Values, class Vectors>
        UTOPIA_INLINE_FUNCTION static bool check(
            const Mat &mat,
            const Values &values,
            const Vectors &vectors)
        {
            auto sum_v = sum(vectors);
            auto sum_e = sum(values);


            bool ok = !device::isnan(sum_v) && !device::isnan(sum_e);

            UTOPIA_DEVICE_ASSERT( !device::isnan(sum_v) );
            UTOPIA_DEVICE_ASSERT( !device::isnan(sum_e) );

            if(ok) {
                //check recomposition
                ok = approxeq(mat, vectors * diag(values) * transpose(vectors), 1e-8);
            }

            if(!ok) {
                disp("--------------------");
                disp(mat);
                disp("--------------------");
                disp(vectors);
                disp("--------------------");
                disp(values);
                disp("--------------------");

                Mat reconstructed = vectors * diag(values) * transpose(vectors);
                disp(reconstructed);
            }

            UTOPIA_DEVICE_ASSERT(ok);
            return ok;
        }

    };

}

#endif //UTOPIA_DEVICE_EIGEN_DECOMPOSITION_HPP
