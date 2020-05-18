#ifndef UTOPIA_ORTHOGONAL_TRANSFORMATION_HPP
#define UTOPIA_ORTHOGONAL_TRANSFORMATION_HPP

#include "moonolith_householder.hpp"
#include "moonolith_vector.hpp"

#include <libmesh/dense_matrix.h>

namespace utopia {

    class OrthogonalTransformation {
    public:
        void build(const std::vector<double> &n, libMesh::DenseMatrix<double> &mat) {
            auto spatial_dim = n.size();

            if (spatial_dim == 2) {
                v2.x = n[0];
                v2.y = n[1];

                build_aux(v2, h2, mat);

            } else if (spatial_dim == 3) {
                v3.x = n[0];
                v3.y = n[1];
                v3.z = n[2];

                build_aux(v3, h3, mat);
            } else {
                assert(false);
            }
        }

    private:
        moonolith::Vector<double, 2> v2;
        moonolith::Vector<double, 3> v3;

        moonolith::HouseholderTransformation<double, 2> h2;
        moonolith::HouseholderTransformation<double, 3> h3;

        template <int Dim>
        static void build_aux(moonolith::Vector<double, Dim> &v,
                              moonolith::HouseholderTransformation<double, Dim> &H,
                              libMesh::DenseMatrix<double> &mat) {
            v[0] -= 1;

            auto len = length(v);

            if (len == 0.) {
                H.identity();
            } else {
                v /= len;
                H.init(v);
            }

            mat.resize(Dim, Dim);

            for (int i = 0; i < Dim; ++i) {
                for (int j = 0; j < Dim; ++j) {
                    mat(i, j) = H(i, j);
                }
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_ORTHOGONAL_TRANSFORMATION_HPP
