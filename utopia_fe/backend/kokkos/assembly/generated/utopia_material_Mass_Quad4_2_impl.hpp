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
                //	- Subexpressions: 10*ADD + 8*DIV + 31*MUL + 4*POW + 15*SUB
                T x0 = pow(1 - x, 2);
                T x1 = pow(1 - y, 2);
                T x2 = (1.0 / 4.0) * x - 0.25;
                T x3 = (1.0 / 4.0) * x + 0.25;
                T x4 = (1.0 / 4.0) * y - 0.25;
                T x5 = (1.0 / 4.0) * y + 0.25;
                T x6 = weight * (-(px[0] * x2 - px[1] * x3 + px[2] * x3 - px[3] * x2) *
                                     (py[0] * x4 - py[1] * x4 + py[2] * x5 - py[3] * x5) +
                                 (px[0] * x4 - px[1] * x4 + px[2] * x5 - px[3] * x5) *
                                     (py[0] * x2 - py[1] * x3 + py[2] * x3 - py[3] * x2));
                T x7 = (1.0 / 16.0) * x6;
                T x8 = x1 * x7;
                T x9 = (1.0 / 4.0) * x6;
                T x10 = (1.0 / 2.0) * x;
                T x11 = (0.5 - x10) * (x10 + 0.5);
                T x12 = x11 * x9;
                T x13 = x1 * x12;
                T x14 = (1.0 / 2.0) * y;
                T x15 = (0.5 - x14) * (x14 + 0.5);
                T x16 = x11 * x15 * x6;
                T x17 = x15 * x9;
                T x18 = x0 * x17;
                T x19 = pow(x + 1, 2);
                T x20 = x17 * x19;
                T x21 = pow(y + 1, 2);
                T x22 = x21 * x7;
                T x23 = x12 * x21;
                H[0] += x0 * x8;
                H[1] += x13;
                H[2] += x16;
                H[3] += x18;
                H[4] += x13;
                H[5] += x19 * x8;
                H[6] += x20;
                H[7] += x16;
                H[8] += x16;
                H[9] += x20;
                H[10] += x19 * x22;
                H[11] += x23;
                H[12] += x18;
                H[13] += x16;
                H[14] += x23;
                H[15] += x0 * x22;
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
                //	- Subexpressions: 11*ADD + 6*DIV + 28*MUL + 13*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x1 * x3;
                T x5 = (1.0 / 4.0) * x - 0.25;
                T x6 = (1.0 / 4.0) * x + 0.25;
                T x7 = (1.0 / 4.0) * y - 0.25;
                T x8 = (1.0 / 4.0) * y + 0.25;
                T x9 = x0 + 0.5;
                T x10 = x3 * x9;
                T x11 = x2 + 0.5;
                T x12 = x11 * x9;
                T x13 = x1 * x11;
                T x14 = weight *
                        (-(px[0] * x5 - px[1] * x6 + px[2] * x6 - px[3] * x5) *
                             (py[0] * x7 - py[1] * x7 + py[2] * x8 - py[3] * x8) +
                         (px[0] * x7 - px[1] * x7 + px[2] * x8 - px[3] * x8) *
                             (py[0] * x5 - py[1] * x6 + py[2] * x6 - py[3] * x5)) *
                        (u[0] * x4 + u[1] * x10 + u[2] * x12 + u[3] * x13);
                Hx[0] += x14 * x4;
                Hx[1] += x10 * x14;
                Hx[2] += x12 * x14;
                Hx[3] += x13 * x14;
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
                //	- Subexpressions: 11*ADD + 6*DIV + 28*MUL + 13*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x1 * x3;
                T x5 = (1.0 / 4.0) * x - 0.25;
                T x6 = (1.0 / 4.0) * x + 0.25;
                T x7 = (1.0 / 4.0) * y - 0.25;
                T x8 = (1.0 / 4.0) * y + 0.25;
                T x9 = x0 + 0.5;
                T x10 = x3 * x9;
                T x11 = x2 + 0.5;
                T x12 = x11 * x9;
                T x13 = x1 * x11;
                T x14 = weight *
                        (-(px[0] * x5 - px[1] * x6 + px[2] * x6 - px[3] * x5) *
                             (py[0] * x7 - py[1] * x7 + py[2] * x8 - py[3] * x8) +
                         (px[0] * x7 - px[1] * x7 + px[2] * x8 - px[3] * x8) *
                             (py[0] * x5 - py[1] * x6 + py[2] * x6 - py[3] * x5)) *
                        (u[0] * x4 + u[1] * x10 + u[2] * x12 + u[3] * x13);
                g[0] += x14 * x4;
                g[1] += x10 * x14;
                g[2] += x12 * x14;
                g[3] += x13 * x14;
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
                //	- Result: 6*ADD + ADDAUGMENTEDASSIGNMENT + 23*MUL + POW
                //	- Subexpressions: 4*ADD + 6*DIV + 4*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x0 + 0.5;
                T x5 = x2 + 0.5;
                T x6 = (1.0 / 4.0) * x - 0.25;
                T x7 = (1.0 / 4.0) * x + 0.25;
                T x8 = (1.0 / 4.0) * y - 0.25;
                T x9 = (1.0 / 4.0) * y + 0.25;
                e += weight *
                     (-(px[0] * x6 - px[1] * x7 + px[2] * x7 - px[3] * x6) *
                          (py[0] * x8 - py[1] * x8 + py[2] * x9 - py[3] * x9) +
                      (px[0] * x8 - px[1] * x8 + px[2] * x9 - px[3] * x9) *
                          (py[0] * x6 - py[1] * x7 + py[2] * x7 - py[3] * x6)) *
                     pow(u[0] * x1 * x3 + u[1] * x3 * x4 + u[2] * x4 * x5 + u[3] * x1 * x5, 2);
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
                //	- Result: 21*ADDAUGMENTEDASSIGNMENT + 9*MUL + POW
                //	- Subexpressions: 13*ADD + 8*DIV + 40*MUL + 4*POW + 15*SUB
                T x0 = pow(1 - x, 2);
                T x1 = pow(1 - y, 2);
                T x2 = (1.0 / 4.0) * x - 0.25;
                T x3 = (1.0 / 4.0) * x + 0.25;
                T x4 = (1.0 / 4.0) * y - 0.25;
                T x5 = (1.0 / 4.0) * y + 0.25;
                T x6 = weight * (-(px[0] * x2 - px[1] * x3 + px[2] * x3 - px[3] * x2) *
                                     (py[0] * x4 - py[1] * x4 + py[2] * x5 - py[3] * x5) +
                                 (px[0] * x4 - px[1] * x4 + px[2] * x5 - px[3] * x5) *
                                     (py[0] * x2 - py[1] * x3 + py[2] * x3 - py[3] * x2));
                T x7 = (1.0 / 16.0) * x6;
                T x8 = x1 * x7;
                T x9 = (1.0 / 2.0) * x;
                T x10 = 0.5 - x9;
                T x11 = x9 + 0.5;
                T x12 = (1.0 / 4.0) * x6;
                T x13 = x10 * x11 * x12;
                T x14 = x1 * x13;
                T x15 = (1.0 / 2.0) * y;
                T x16 = 0.5 - x15;
                T x17 = x10 * x16;
                T x18 = x15 + 0.5;
                T x19 = x11 * x18;
                T x20 = x17 * x19 * x6;
                T x21 = x12 * x16 * x18;
                T x22 = x0 * x21;
                T x23 = pow(x + 1, 2);
                T x24 = x21 * x23;
                T x25 = pow(y + 1, 2);
                T x26 = x25 * x7;
                T x27 = x13 * x25;
                T x28 = x11 * x16;
                T x29 = x10 * x18;
                T x30 = u[0] * x17 + u[1] * x28 + u[2] * x19 + u[3] * x29;
                T x31 = x30 * x6;
                H[0] += x0 * x8;
                H[1] += x14;
                H[2] += x20;
                H[3] += x22;
                H[4] += x14;
                H[5] += x23 * x8;
                H[6] += x24;
                H[7] += x20;
                H[8] += x20;
                H[9] += x24;
                H[10] += x23 * x26;
                H[11] += x27;
                H[12] += x22;
                H[13] += x20;
                H[14] += x27;
                H[15] += x0 * x26;
                g[0] += x17 * x31;
                g[1] += x28 * x31;
                g[2] += x19 * x31;
                g[3] += x29 * x31;
                e += pow(x30, 2) * x6;
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
