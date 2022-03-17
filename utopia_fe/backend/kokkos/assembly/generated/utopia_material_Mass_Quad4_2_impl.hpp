#ifndef UTOPIA_TPL_MATERIAL_Mass_Quad4_2_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_Quad4_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Quad4_2.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Mass for symmetric element pair trial=test=Quad4
         */
        template <typename T, typename GeoT>
        class Mass<Quad4<T, GeoT>> {
        public:
            using ElemT = Quad4<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Mass<Quad4>"; }

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
                //	- Result: 16*ADDAUGMENTEDASSIGNMENT + 4*MUL
                //	- Subexpressions: 15*ADD + 4*DIV + 66*MUL + 4*POW + 16*SUB
                T x0 = pow(1 - x, 2);
                T x1 = pow(1 - y, 2);
                T x2 = 0.125 * px[0];
                T x3 = 0.125 * px[1];
                T x4 = 0.125 * px[2];
                T x5 = py[1] * x4;
                T x6 = 0.125 * px[3];
                T x7 = py[2] * x6;
                T x8 = py[3] * x3;
                T x9 = weight * (0.125 * px[0] * py[1] + 0.125 * px[0] * py[2] * y + 0.125 * px[0] * py[3] * x +
                                 0.125 * px[1] * py[0] * y + 0.125 * px[1] * py[2] * x + 0.125 * px[1] * py[2] +
                                 0.125 * px[2] * py[0] * x + 0.125 * px[2] * py[3] * y + 0.125 * px[2] * py[3] +
                                 0.125 * px[3] * py[0] + 0.125 * px[3] * py[1] * x + 0.125 * px[3] * py[1] * y -
                                 py[0] * x * x6 - py[0] * x3 - py[0] * x4 * y - py[1] * x2 * y - py[2] * x * x2 -
                                 py[3] * x2 - x * x5 - x * x8 - x5 - x7 * y - x7 - x8 * y);
                T x10 = (1.0 / 16.0) * x9;
                T x11 = x1 * x10;
                T x12 = (1.0 / 4.0) * x9;
                T x13 = (1.0 / 2.0) * x;
                T x14 = (0.5 - x13) * (x13 + 0.5);
                T x15 = x12 * x14;
                T x16 = x1 * x15;
                T x17 = (1.0 / 2.0) * y;
                T x18 = (0.5 - x17) * (x17 + 0.5);
                T x19 = x14 * x18 * x9;
                T x20 = x12 * x18;
                T x21 = x0 * x20;
                T x22 = pow(x + 1, 2);
                T x23 = x20 * x22;
                T x24 = pow(y + 1, 2);
                T x25 = x10 * x24;
                T x26 = x15 * x24;
                H[0] += x0 * x11;
                H[1] += x16;
                H[2] += x19;
                H[3] += x21;
                H[4] += x16;
                H[5] += x11 * x22;
                H[6] += x23;
                H[7] += x19;
                H[8] += x19;
                H[9] += x23;
                H[10] += x22 * x25;
                H[11] += x26;
                H[12] += x21;
                H[13] += x19;
                H[14] += x26;
                H[15] += x0 * x25;
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
                //	- Result: 4*ADDAUGMENTEDASSIGNMENT + 4*MUL
                //	- Subexpressions: 16*ADD + 2*DIV + 63*MUL + 14*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x1 * x3;
                T x5 = x0 + 0.5;
                T x6 = x3 * x5;
                T x7 = x2 + 0.5;
                T x8 = x5 * x7;
                T x9 = x1 * x7;
                T x10 = 0.125 * px[0];
                T x11 = 0.125 * px[1];
                T x12 = 0.125 * px[2];
                T x13 = py[1] * x12;
                T x14 = 0.125 * px[3];
                T x15 = py[2] * x14;
                T x16 = py[3] * x11;
                T x17 = weight * (u[0] * x4 + u[1] * x6 + u[2] * x8 + u[3] * x9) *
                        (0.125 * px[0] * py[1] + 0.125 * px[0] * py[2] * y + 0.125 * px[0] * py[3] * x +
                         0.125 * px[1] * py[0] * y + 0.125 * px[1] * py[2] * x + 0.125 * px[1] * py[2] +
                         0.125 * px[2] * py[0] * x + 0.125 * px[2] * py[3] * y + 0.125 * px[2] * py[3] +
                         0.125 * px[3] * py[0] + 0.125 * px[3] * py[1] * x + 0.125 * px[3] * py[1] * y -
                         py[0] * x * x14 - py[0] * x11 - py[0] * x12 * y - py[1] * x10 * y - py[2] * x * x10 -
                         py[3] * x10 - x * x13 - x * x16 - x13 - x15 * y - x15 - x16 * y);
                Hx[0] += x17 * x4;
                Hx[1] += x17 * x6;
                Hx[2] += x17 * x8;
                Hx[3] += x17 * x9;
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
                //	- Result: 4*ADDAUGMENTEDASSIGNMENT + 4*MUL
                //	- Subexpressions: 16*ADD + 2*DIV + 63*MUL + 14*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x1 * x3;
                T x5 = x0 + 0.5;
                T x6 = x3 * x5;
                T x7 = x2 + 0.5;
                T x8 = x5 * x7;
                T x9 = x1 * x7;
                T x10 = 0.125 * px[0];
                T x11 = 0.125 * px[1];
                T x12 = 0.125 * px[2];
                T x13 = py[1] * x12;
                T x14 = 0.125 * px[3];
                T x15 = py[2] * x14;
                T x16 = py[3] * x11;
                T x17 = weight * (u[0] * x4 + u[1] * x6 + u[2] * x8 + u[3] * x9) *
                        (0.125 * px[0] * py[1] + 0.125 * px[0] * py[2] * y + 0.125 * px[0] * py[3] * x +
                         0.125 * px[1] * py[0] * y + 0.125 * px[1] * py[2] * x + 0.125 * px[1] * py[2] +
                         0.125 * px[2] * py[0] * x + 0.125 * px[2] * py[3] * y + 0.125 * px[2] * py[3] +
                         0.125 * px[3] * py[0] + 0.125 * px[3] * py[1] * x + 0.125 * px[3] * py[1] * y -
                         py[0] * x * x14 - py[0] * x11 - py[0] * x12 * y - py[1] * x10 * y - py[2] * x * x10 -
                         py[3] * x10 - x * x13 - x * x16 - x13 - x15 * y - x15 - x16 * y);
                g[0] += x17 * x4;
                g[1] += x17 * x6;
                g[2] += x17 * x8;
                g[3] += x17 * x9;
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
                // TODO
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
                //	- Result: 20*ADDAUGMENTEDASSIGNMENT + 8*MUL
                //	- Subexpressions: 18*ADD + 4*DIV + 75*MUL + 4*POW + 16*SUB
                T x0 = pow(1 - x, 2);
                T x1 = pow(1 - y, 2);
                T x2 = 0.125 * px[0];
                T x3 = 0.125 * px[1];
                T x4 = 0.125 * px[2];
                T x5 = py[1] * x4;
                T x6 = 0.125 * px[3];
                T x7 = py[2] * x6;
                T x8 = py[3] * x3;
                T x9 = weight * (0.125 * px[0] * py[1] + 0.125 * px[0] * py[2] * y + 0.125 * px[0] * py[3] * x +
                                 0.125 * px[1] * py[0] * y + 0.125 * px[1] * py[2] * x + 0.125 * px[1] * py[2] +
                                 0.125 * px[2] * py[0] * x + 0.125 * px[2] * py[3] * y + 0.125 * px[2] * py[3] +
                                 0.125 * px[3] * py[0] + 0.125 * px[3] * py[1] * x + 0.125 * px[3] * py[1] * y -
                                 py[0] * x * x6 - py[0] * x3 - py[0] * x4 * y - py[1] * x2 * y - py[2] * x * x2 -
                                 py[3] * x2 - x * x5 - x * x8 - x5 - x7 * y - x7 - x8 * y);
                T x10 = (1.0 / 16.0) * x9;
                T x11 = x1 * x10;
                T x12 = (1.0 / 2.0) * x;
                T x13 = 0.5 - x12;
                T x14 = x12 + 0.5;
                T x15 = (1.0 / 4.0) * x9;
                T x16 = x13 * x14 * x15;
                T x17 = x1 * x16;
                T x18 = (1.0 / 2.0) * y;
                T x19 = 0.5 - x18;
                T x20 = x13 * x19;
                T x21 = x18 + 0.5;
                T x22 = x14 * x21;
                T x23 = x20 * x22 * x9;
                T x24 = x15 * x19 * x21;
                T x25 = x0 * x24;
                T x26 = pow(x + 1, 2);
                T x27 = x24 * x26;
                T x28 = pow(y + 1, 2);
                T x29 = x10 * x28;
                T x30 = x16 * x28;
                T x31 = x14 * x19;
                T x32 = x13 * x21;
                T x33 = x9 * (u[0] * x20 + u[1] * x31 + u[2] * x22 + u[3] * x32);
                H[0] += x0 * x11;
                H[1] += x17;
                H[2] += x23;
                H[3] += x25;
                H[4] += x17;
                H[5] += x11 * x26;
                H[6] += x27;
                H[7] += x23;
                H[8] += x23;
                H[9] += x27;
                H[10] += x26 * x29;
                H[11] += x30;
                H[12] += x25;
                H[13] += x23;
                H[14] += x30;
                H[15] += x0 * x29;
                g[0] += x20 * x33;
                g[1] += x31 * x33;
                g[2] += x22 * x33;
                g[3] += x32 * x33;
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using MassQuad4 = utopia::kokkos::
            AutoKernel<FE, utopia::kernels::Mass<utopia::kernels::Quad4<typename FE::Scalar, typename FE::Scalar>>, 2>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_Mass_Quad4_2_IMPL_hpp
