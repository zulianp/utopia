#ifndef UTOPIA_ELEMENT_WISE_PSEUDO_INVERSE_HPP
#define UTOPIA_ELEMENT_WISE_PSEUDO_INVERSE_HPP

#include "utopia_Base.hpp"

#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"

#include <cmath>
#include <limits>

namespace utopia {

    template <class Vector, int Backend = Traits<Vector>::Backend>
    class ElementWisePseudoInverse {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        static void apply(const Vector &in, Vector &out, const Scalar tol = std::numeric_limits<Scalar>::epsilon()) {
            if (in.is_alias(out)) {
                auto out_view = local_view_device(out);

                parallel_for(local_range_device(out), UTOPIA_LAMBDA(const SizeType &i) {
                    auto val = out_view.get(i);
                    if (device::abs(val) > tol) {
                        val = 1. / val;
                    } else {
                        val = 0.0;
                    }

                    out_view.set(i, val);
                });

            } else {
                auto in_view = const_local_view_device(in);

                if (empty(out) || out.size() != in.size()) {
                    out.zeros(layout(in));
                }

                auto out_view = local_view_device(out);

                parallel_for(local_range_device(out), UTOPIA_LAMBDA(const SizeType &i) {
                    auto val = in_view.get(i);
                    if (device::abs(val) > tol) {
                        val = 1. / val;
                    } else {
                        val = 0.0;
                    }

                    out_view.set(i, val);
                });
            }
        }
    };

    template <class Vector>
    void e_pseudo_inv(const Vector &in,
                      Vector &out,
                      const UTOPIA_SCALAR(Vector) tol = std::numeric_limits<UTOPIA_SCALAR(Vector)>::epsilon()) {
        ElementWisePseudoInverse<Vector>::apply(in, out, tol);
    }
}  // namespace utopia

#endif  // UTOPIA_ELEMENT_WISE_PSEUDO_INVERSE_HPP
