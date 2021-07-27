#include "utopia_moonolith_HouseholderReflection.hpp"

#include "moonolith_householder.hpp"
#include "utopia.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Dim>
    void HouseholderReflectionForContact<Matrix, Vector, Dim>::build(const Vector &is_contact,
                                                                     const Vector &normal,
                                                                     Matrix &trafo) {
        using Scalar = typename Traits<Vector>::Scalar;

        trafo.sparse(square_matrix_layout(layout(normal)), Dim, 0);

        ::moonolith::HouseholderTransformation<Scalar, Dim> H;
        auto r = range(normal);

        Read<Vector> r_normal(normal), r_ic(is_contact);
        Write<Matrix> w_ot(trafo, utopia::LOCAL);

        ::moonolith::Vector<Scalar, Dim> n;
        for (auto i = r.begin(); i < r.end(); i += Dim) {
            if (is_contact.get(i) == 0.0) {
                H.identity();
            } else {
                for (int d = 0; d < Dim; ++d) {
                    n[d] = normal.get(i + d);
                }

                auto len = length(n);
                assert(len > 0.0);

                n /= len;
                n.x -= 1.0;

                len = length(n);

                if (approxeq(len, 0.0, 1e-15)) {
                    H.identity();
                } else {
                    n /= len;
                    H.init(n);
                }
            }

            assert(approxeq(std::abs(measure(H)), 1.0, 1e-8));

            for (int d1 = 0; d1 < Dim; ++d1) {
                for (int d2 = 0; d2 < Dim; ++d2) {
                    trafo.set(i + d1, i + d2, H(d1, d2));
                }
            }
        }
    }
}  // namespace utopia
