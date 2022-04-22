#ifndef UTOPIA_TPL_MATERIAL_Mass_Tri3_2_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_Tri3_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Tri3_2.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Mass for symmetric element pair trial=test=Tri3
         */
        template <typename T, typename GeoT>
        class Mass<Tri3<T, GeoT>> {
        public:
            using ElemT = Tri3<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Mass<Tri3>"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    // TODO
                }

                // TODO
            };

            Mass(const Params &params = Params()) {
                // TODO
            }

            UTOPIA_FUNCTION void hessian(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 9*ADDAUGMENTEDASSIGNMENT + 3*MUL + 3*POW
                //	- Subexpressions: 8*MUL + 7*SUB
                T x0 = -x - y + 1;
                T x1 = weight * ((-px[0] + px[1]) * (-py[0] + py[2]) - (-px[0] + px[2]) * (-py[0] + py[1]));
                T x2 = x0 * x1;
                T x3 = x * x2;
                T x4 = x2 * y;
                T x5 = x * x1 * y;
                H[0] += pow(x0, 2) * x1;
                H[1] += x3;
                H[2] += x4;
                H[3] += x3;
                H[4] += pow(x, 2) * x1;
                H[5] += x5;
                H[6] += x4;
                H[7] += x5;
                H[8] += x1 * pow(y, 2);
            }

            UTOPIA_FUNCTION void apply(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T *UTOPIA_RESTRICT Hx) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 3*ADDAUGMENTEDASSIGNMENT + 3*MUL
                //	- Subexpressions: 2*ADD + 7*MUL + 7*SUB
                T x0 = -x - y + 1;
                T x1 = weight * ((-px[0] + px[1]) * (-py[0] + py[2]) - (-px[0] + px[2]) * (-py[0] + py[1])) *
                       (u[0] * x0 + u[1] * x + u[2] * y);
                Hx[0] += x0 * x1;
                Hx[1] += x * x1;
                Hx[2] += x1 * y;
            }

            UTOPIA_FUNCTION void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T *UTOPIA_RESTRICT g) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 3*ADDAUGMENTEDASSIGNMENT + 3*MUL
                //	- Subexpressions: 2*ADD + 7*MUL + 7*SUB
                T x0 = -x - y + 1;
                T x1 = weight * ((-px[0] + px[1]) * (-py[0] + py[2]) - (-px[0] + px[2]) * (-py[0] + py[1])) *
                       (u[0] * x0 + u[1] * x + u[2] * y);
                g[0] += x0 * x1;
                g[1] += x * x1;
                g[2] += x1 * y;
            }

            UTOPIA_FUNCTION void value(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T &e) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 7*ADD + ADDAUGMENTEDASSIGNMENT + 12*MUL + POW
                //	- Subexpressions: 0
                e += weight * ((-px[0] + px[1]) * (-py[0] + py[2]) - (-px[0] + px[2]) * (-py[0] + py[1])) *
                     pow(u[0] * (-x - y + 1) + u[1] * x + u[2] * y, 2);
            }

            UTOPIA_FUNCTION void eval(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T &e,
                T *UTOPIA_RESTRICT g,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 13*ADDAUGMENTEDASSIGNMENT + 7*MUL + 4*POW
                //	- Subexpressions: 2*ADD + 11*MUL + 7*SUB
                T x0 = -x - y + 1;
                T x1 = weight * ((-px[0] + px[1]) * (-py[0] + py[2]) - (-px[0] + px[2]) * (-py[0] + py[1]));
                T x2 = x0 * x1;
                T x3 = x * x2;
                T x4 = x2 * y;
                T x5 = x * x1;
                T x6 = x5 * y;
                T x7 = u[0] * x0 + u[1] * x + u[2] * y;
                H[0] += pow(x0, 2) * x1;
                H[1] += x3;
                H[2] += x4;
                H[3] += x3;
                H[4] += pow(x, 2) * x1;
                H[5] += x6;
                H[6] += x4;
                H[7] += x6;
                H[8] += x1 * pow(y, 2);
                g[0] += x2 * x7;
                g[1] += x5 * x7;
                g[2] += x1 * x7 * y;
                e += x1 * pow(x7, 2);
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using MassTri3 = utopia::kokkos::
            AutoKernel<FE, utopia::kernels::Mass<utopia::kernels::Tri3<typename FE::Scalar, typename FE::Scalar>>, 2>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_Mass_Tri3_2_IMPL_hpp
