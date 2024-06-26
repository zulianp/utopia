#ifndef UTOPIA_TPL_FE_Quad4_2_IMPL_hpp
#define UTOPIA_TPL_FE_Quad4_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

// #include "utopia_fe_Quad4.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Quad4 for dimension 2
         */
        template <typename T, typename GeoT>
        class Quad4 {
        public:
            static constexpr int Dim = 2;
            static constexpr int NNodes = 4;
            static constexpr int Order = 1;

            using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Quad4"; }

            UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

            UTOPIA_INLINE_FUNCTION static constexpr int n_nodes() { return NNodes; }

            UTOPIA_INLINE_FUNCTION static constexpr int order() { return Order; }

            UTOPIA_FUNCTION static Result measure(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Input quadrature point
                const T x,
                const T y) {
                T measure_value = 0;
                // FLOATING POINT OPS!
                //	- Result: 5*ADD + ASSIGNMENT + 18*MUL
                //	- Subexpressions: 2*ADD + 4*DIV + 2*SUB
                T x0 = (1.0 / 4.0) * x - 0.25;
                T x1 = (1.0 / 4.0) * x + 0.25;
                T x2 = (1.0 / 4.0) * y - 0.25;
                T x3 = (1.0 / 4.0) * y + 0.25;
                measure_value = -(px[0] * x0 - px[1] * x1 + px[2] * x1 - px[3] * x0) *
                                    (py[0] * x2 - py[1] * x2 + py[2] * x3 - py[3] * x3) +
                                (px[0] * x2 - px[1] * x2 + px[2] * x3 - px[3] * x3) *
                                    (py[0] * x0 - py[1] * x1 + py[2] * x1 - py[3] * x0);
                return measure_value;
            }

            UTOPIA_FUNCTION static void jacobian(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Input quadrature point
                const T x,
                const T y,
                GeoT *UTOPIA_RESTRICT J) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 4*ADD + 4*ASSIGNMENT + 16*MUL
                //	- Subexpressions: 2*ADD + 4*DIV + 2*SUB
                T x0 = (1.0 / 4.0) * y - 0.25;
                T x1 = (1.0 / 4.0) * y + 0.25;
                T x2 = (1.0 / 4.0) * x - 0.25;
                T x3 = (1.0 / 4.0) * x + 0.25;
                J[0] = px[0] * x0 - px[1] * x0 + px[2] * x1 - px[3] * x1;
                J[1] = px[0] * x2 - px[1] * x3 + px[2] * x3 - px[3] * x2;
                J[2] = py[0] * x0 - py[1] * x0 + py[2] * x1 - py[3] * x1;
                J[3] = py[0] * x2 - py[1] * x3 + py[2] * x3 - py[3] * x2;
            }

            UTOPIA_FUNCTION static void jacobian_inverse(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Input quadrature point
                const T x,
                const T y,
                GeoT *UTOPIA_RESTRICT J_inv) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 4*ASSIGNMENT + 4*MUL
                //	- Subexpressions: 6*ADD + DIV + 18*MUL + 11*SUB
                T x0 = x - 1.0;
                T x1 = x + 1.0;
                T x2 = py[0] * x0 - py[1] * x1 + py[2] * x1 - py[3] * x0;
                T x3 = px[0] * x0 - px[1] * x1 + px[2] * x1 - px[3] * x0;
                T x4 = y - 1.0;
                T x5 = y + 1.0;
                T x6 = py[0] * x4 - py[1] * x4 + py[2] * x5 - py[3] * x5;
                T x7 = px[0] * x4 - px[1] * x4 + px[2] * x5 - px[3] * x5;
                T x8 = 4 / (-x2 * x7 + x3 * x6);
                J_inv[0] = -x2 * x8;
                J_inv[1] = x3 * x8;
                J_inv[2] = x6 * x8;
                J_inv[3] = -x7 * x8;
            }

            UTOPIA_FUNCTION static void transform(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Input quadrature point
                const T x,
                const T y,
                GeoT &tx,
                GeoT &ty) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 2*ADD + 2*ASSIGNMENT + 8*MUL
                //	- Subexpressions: 2*ADD + 2*DIV + 4*MUL + 2*SUB
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
                tx = px[0] * x4 + px[1] * x6 + px[2] * x8 + px[3] * x9;
                ty = py[0] * x4 + py[1] * x6 + py[2] * x8 + py[3] * x9;
            }

            UTOPIA_FUNCTION static void inverse_transform(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Input quadrature point
                const T tx,
                const T ty,
                GeoT &x,
                GeoT &y) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 4*ADD + 2*ASSIGNMENT + 23*MUL + 2*POW
                //	- Subexpressions: 49*ADD + DIV + 143*MUL + 11*POW + 41*SUB
                T x0 = px[0] * py[2];
                T x1 = px[3] * py[0];
                T x2 = px[0] * py[3];
                T x3 = px[2] * py[0];
                T x4 = px[1] * py[3];
                T x5 = px[3] * py[1];
                T x6 = x4 - x5;
                T x7 = px[2] * py[1];
                T x8 = px[1] * py[2];
                T x9 = x7 - x8;
                T x10 = 0.5 * x0;
                T x11 = px[1] * ty;
                T x12 = px[0] * ty;
                T x13 = 0.5 * x12;
                T x14 = px[2] * ty;
                T x15 = px[3] * ty;
                T x16 = py[0] * tx;
                T x17 = 0.5 * x16;
                T x18 = py[1] * tx;
                T x19 = 0.5 * x18;
                T x20 = py[3] * tx;
                T x21 = 0.5 * x20;
                T x22 = 0.5 * tx;
                T x23 = ty * x22;
                T x24 = 0.5 * x4;
                T x25 = 0.5 * x11;
                T x26 = py[2] * tx;
                T x27 = 0.5 * x26;
                T x28 = 0.5 * x1;
                T x29 = pow(ty, 2);
                T x30 = 0.5 * x29;
                T x31 = px[0] * x30;
                T x32 = pow(py[2], 2);
                T x33 = pow(px[0], 2);
                T x34 = 0.5 * ty;
                T x35 = px[1] * x30;
                T x36 = pow(py[3], 2);
                T x37 = pow(px[1], 2);
                T x38 = pow(py[0], 2);
                T x39 = pow(px[2], 2);
                T x40 = pow(py[1], 2);
                T x41 = pow(px[3], 2);
                T x42 = pow(tx, 2);
                T x43 = 0.5 * x42;
                T x44 = py[0] * x43;
                T x45 = py[1] * x43;
                T x46 = 0.25 * x29;
                T x47 = 0.25 * x42;
                T x48 =
                    2.0 *
                    sqrt(-px[0] * x22 * x32 - px[1] * x22 * x36 - px[1] * x31 - px[2] * px[3] * x30 -
                         px[2] * x22 * x38 + px[2] * x31 - px[2] * x35 - px[3] * x22 * x40 - px[3] * x31 + px[3] * x35 -
                         py[0] * x34 * x39 - py[1] * x34 * x41 - py[1] * x44 - py[2] * py[3] * x43 - py[2] * x33 * x34 +
                         py[2] * x44 - py[2] * x45 - py[3] * x34 * x37 - py[3] * x44 + py[3] * x45 + x0 * x17 +
                         x0 * x19 + x0 * x21 - x0 * x23 - x1 * x11 + x1 * x19 + x1 * x23 - x1 * x26 + x1 * x8 +
                         x10 * x11 + x10 * x14 + x10 * x15 - x10 * x3 - x10 * x4 - x10 * x5 + x11 * x17 - x11 * x19 -
                         x12 * x17 + x12 * x19 - x12 * x7 + x13 * x3 + x13 * x4 + x13 * x5 - x14 * x2 + x14 * x21 +
                         x14 * x24 - x14 * x27 + x14 * x28 - x15 * x21 + x15 * x24 + x15 * x27 + 0.5 * x15 * x7 -
                         x15 * x8 - x16 * x8 + x17 * x4 + x17 * x7 - x18 * x2 + x19 * x4 + x2 * x23 + x2 * x7 -
                         x20 * x7 + x21 * x3 + x21 * x5 - x23 * x3 - x23 * x4 - x23 * x5 + x23 * x7 + x23 * x8 -
                         x24 * x3 - x24 * x5 + x25 * x3 + x25 * x5 + x27 * x3 + x27 * x4 + x27 * x5 - x28 * x7 +
                         0.25 * x32 * x33 + x32 * x47 + x33 * x46 + 0.25 * x36 * x37 + x36 * x47 + x37 * x46 +
                         0.25 * x38 * x39 + x38 * x47 + x39 * x46 + 0.25 * x40 * x41 + x40 * x47 + x41 * x46);
                T x49 = px[1] * py[0];
                T x50 = -px[2] * py[3] + px[3] * py[2];
                x = (px[3] * py[0] - x11 + x12 + x14 - x15 - x16 + x18 - x2 + x20 - x26 - x48 - x9) /
                    (x0 + x1 - x2 - x3 + x6 + x9);
                y = (px[0] * py[1] + x11 - x12 - x14 + x15 + x16 - x18 - x20 + x26 - x48 - x49 - x50) /
                    (px[0] * py[1] - x0 + x3 - x49 + x50 + x6);
            }

            UTOPIA_FUNCTION static void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Input quadrature point
                const T x,
                const T y,
                // Output
                Result *UTOPIA_RESTRICT gx,
                Result *UTOPIA_RESTRICT gy) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 8*ADD + 8*ASSIGNMENT + 16*MUL
                //	- Subexpressions: 6*ADD + DIV + 26*MUL + 4*NEG + 11*SUB
                T x0 = x - 1.0;
                T x1 = x + 1.0;
                T x2 = px[0] * x0 - px[1] * x1 + px[2] * x1 - px[3] * x0;
                T x3 = x0 * x2;
                T x4 = y - 1.0;
                T x5 = py[0] * x0 - py[1] * x1 + py[2] * x1 - py[3] * x0;
                T x6 = -x4 * x5;
                T x7 = y + 1.0;
                T x8 = py[0] * x4 - py[1] * x4 + py[2] * x7 - py[3] * x7;
                T x9 = px[0] * x4 - px[1] * x4 + px[2] * x7 - px[3] * x7;
                T x10 = 1.0 / (x2 * x8 - x5 * x9);
                T x11 = x0 * x9;
                T x12 = -x4 * x8;
                T x13 = x1 * x2;
                T x14 = x1 * x9;
                T x15 = -x5 * x7;
                T x16 = -x7 * x8;
                gx[0] = x10 * (x3 + x6);
                gy[0] = x10 * (-x11 - x12);
                gx[1] = x10 * (-x13 - x6);
                gy[1] = x10 * (x12 + x14);
                gx[2] = x10 * (x13 + x15);
                gy[2] = x10 * (-x14 - x16);
                gx[3] = x10 * (-x15 - x3);
                gy[3] = x10 * (x11 + x16);
            }

            UTOPIA_FUNCTION static void value(const T x, const T y, Result *UTOPIA_RESTRICT f) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 4*ASSIGNMENT + 4*MUL
                //	- Subexpressions: 2*ADD + 2*DIV + 2*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x0 + 0.5;
                T x5 = x2 + 0.5;
                f[0] = x1 * x3;
                f[1] = x3 * x4;
                f[2] = x4 * x5;
                f[3] = x1 * x5;
            }

            UTOPIA_FUNCTION static void eval(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Input quadrature point
                const T x,
                const T y,
                // Output
                Result *UTOPIA_RESTRICT f,
                Result *UTOPIA_RESTRICT gx,
                Result *UTOPIA_RESTRICT gy,
                T &measure_value) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 13*ADD + 13*ASSIGNMENT + 38*MUL
                //	- Subexpressions: 8*ADD + 3*DIV + 26*MUL + 4*NEG + 13*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x0 + 0.5;
                T x5 = x2 + 0.5;
                T x6 = x - 1.0;
                T x7 = px[0] * x6;
                T x8 = x + 1.0;
                T x9 = px[1] * x8;
                T x10 = px[2] * x8;
                T x11 = px[3] * x6;
                T x12 = y - 1.0;
                T x13 = py[0] * x12;
                T x14 = py[1] * x12;
                T x15 = y + 1.0;
                T x16 = py[2] * x15;
                T x17 = py[3] * x15;
                T x18 = px[0] * x12;
                T x19 = px[1] * x12;
                T x20 = px[2] * x15;
                T x21 = px[3] * x15;
                T x22 = py[0] * x6;
                T x23 = py[1] * x8;
                T x24 = py[2] * x8;
                T x25 = py[3] * x6;
                T x26 = x10 - x11 + x7 - x9;
                T x27 = x26 * x6;
                T x28 = x22 - x23 + x24 - x25;
                T x29 = -x12 * x28;
                T x30 = x13 - x14 + x16 - x17;
                T x31 = x18 - x19 + x20 - x21;
                T x32 = 1.0 / (x26 * x30 - x28 * x31);
                T x33 = x31 * x6;
                T x34 = -x12 * x30;
                T x35 = x26 * x8;
                T x36 = x31 * x8;
                T x37 = -x15 * x28;
                T x38 = -x15 * x30;
                f[0] = x1 * x3;
                f[1] = x3 * x4;
                f[2] = x4 * x5;
                f[3] = x1 * x5;
                measure_value = -((1.0 / 4.0) * x10 - 1.0 / 4.0 * x11 + (1.0 / 4.0) * x7 - 1.0 / 4.0 * x9) *
                                    ((1.0 / 4.0) * x13 - 1.0 / 4.0 * x14 + (1.0 / 4.0) * x16 - 1.0 / 4.0 * x17) +
                                ((1.0 / 4.0) * x18 - 1.0 / 4.0 * x19 + (1.0 / 4.0) * x20 - 1.0 / 4.0 * x21) *
                                    ((1.0 / 4.0) * x22 - 1.0 / 4.0 * x23 + (1.0 / 4.0) * x24 - 1.0 / 4.0 * x25);
                gx[0] = x32 * (x27 + x29);
                gy[0] = x32 * (-x33 - x34);
                gx[1] = x32 * (-x29 - x35);
                gy[1] = x32 * (x34 + x36);
                gx[2] = x32 * (x35 + x37);
                gy[2] = x32 * (-x36 - x38);
                gx[3] = x32 * (-x27 - x37);
                gy[3] = x32 * (x33 + x38);
            }
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_FE_Quad4_2_IMPL_hpp
