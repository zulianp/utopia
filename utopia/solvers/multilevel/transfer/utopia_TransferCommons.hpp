#ifndef UTOPIA_TRANSFER_COMMONS_HPP
#define UTOPIA_TRANSFER_COMMONS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

#include <limits>

namespace utopia {

    template <class Matrix, class Vector>
    class BooleanRestrictOR {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        static bool apply(const Matrix &R, const Vector &x, Vector &x_new) {
            static const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

            Matrix R_boolean = R;

            R_boolean.transform_values(UTOPIA_LAMBDA(const Scalar &value)->Scalar {
                if (device::abs(value) > off_diag_tol) {
                    return 1.;
                } else {
                    return 0.;
                }
            });

            x_new = R_boolean * x;

            x_new.transform_values(UTOPIA_LAMBDA(const Scalar &value)->Scalar {
                if (value > 1.) {
                    return 1.0;
                } else {
                    return 0.0;
                }
            });

            return true;
        }
    };

    template <class Matrix, class Vector>
    bool boolean_restrict_or(const Matrix &R, const Vector &x, Vector &x_new) {
        return BooleanRestrictOR<Matrix, Vector>::apply(R, x, x_new);
    }

}  // namespace utopia

#endif  // UTOPIA_TRANSFER_COMMONS_HPP
