#ifndef UTOPIA_ELEMENT_WISE_PSEUDO_INVERSE_HPP
#define UTOPIA_ELEMENT_WISE_PSEUDO_INVERSE_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

#include <limits>

namespace utopia {

    template<class Vector, int Backend = Traits<Vector>::Backend>
    class ElementWisePseudoInverse {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        static void apply(const Vector &in, Vector &out, const Scalar tol = std::numeric_limits<Scalar>::epsilon())
        {
            each_transform(in, out, [tol](const SizeType i, const Scalar value) -> Scalar {\
                UTOPIA_UNUSED(i);

                if(std::abs(value) > tol) {
                    return 1./value;
                } else {
                    return 0.0;
                }
            });
        }
    };

    template<class Vector>
    void e_pseudo_inv(const Vector &in,  Vector &out, const UTOPIA_SCALAR(Vector) tol = std::numeric_limits<UTOPIA_SCALAR(Vector)>::epsilon())
    {
        ElementWisePseudoInverse<Vector>::apply(in, out, tol);
    }
}

#endif //UTOPIA_ELEMENT_WISE_PSEUDO_INVERSE_HPP
