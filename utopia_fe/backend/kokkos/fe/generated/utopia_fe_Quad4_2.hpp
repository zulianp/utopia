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
        template <typename T, typename GeoT = T>
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

            UTOPIA_FUNCTION static constexpr Result measure(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Input quadrature point
                const T x,
                const T y) {
                T measure_value;
                // FLOATING POINT OPS!
                //	- Result: ADD + ASSIGNMENT + 16*MUL
                //	- Subexpressions: 0
                measure_value = -px[0] * py[1] * y + px[0] * py[1] - px[0] * py[2] * x + px[0] * py[2] * y +
                                px[0] * py[3] * x - px[0] * py[3] * y + px[1] * py[0] * y - px[1] * py[0] +
                                px[1] * py[2] * x - px[1] * py[3] * x + px[2] * py[0] * x - px[2] * py[0] * y -
                                px[2] * py[1] * x - px[3] * py[0] * x + px[3] * py[0] * y + px[3] * py[1] * x;
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
                //	- Subexpressions: 2*SUB
                T x0 = y - 1;
                T x1 = x - 1;
                J[0] = px[0] * x0 - px[1] * x0 + px[2] * y - px[3] * y;
                J[1] = px[0] * x1 - px[1] * x + px[2] * x - px[3] * x;
                J[2] = py[0] * x0 - py[1] * x0 + py[2] * y - py[3] * y;
                J[3] = py[0] * x1 - py[1] * x + py[2] * x - py[3] * x;
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
                //	- Result: 4*ADD + 4*ASSIGNMENT + 19*MUL
                //	- Subexpressions: 7*ADD + DIV + 23*MUL + 8*SUB
                T x0 = py[1] * x;
                T x1 = py[3] * x;
                T x2 = py[0] * x;
                T x3 = py[2] * x;
                T x4 = px[1] * py[0];
                T x5 = px[0] * py[1];
                T x6 = py[3] * y;
                T x7 = py[0] * y;
                T x8 = py[2] * y;
                T x9 =
                    1.0 / (-px[0] * x1 + px[0] * x3 + px[0] * x6 - px[0] * x8 + px[1] * x1 - px[1] * x3 + px[2] * x0 -
                           px[2] * x2 + px[2] * x7 - px[3] * x0 + px[3] * x2 - px[3] * x7 - x4 * y + x4 + x5 * y - x5);
                J_inv[0] = x9 * (py[0] + x0 + x1 - x2 - x3);
                J_inv[1] = x9 * (px[0] * x - px[0] - px[1] * x + px[2] * x - px[3] * x);
                J_inv[2] = x9 * (-py[0] - py[1] * y + py[1] - x6 + x7 + x8);
                J_inv[3] = x9 * (-px[0] * y + px[0] + px[1] * y - px[1] - px[2] * y + px[3] * y);
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
                //	- Subexpressions: 3*MUL + 3*SUB
                T x0 = x * y;
                T x1 = 1 - x0;
                T x2 = 1 - y;
                T x3 = x * x2;
                T x4 = x2 * (1 - x);
                tx = px[0] * x4 + px[1] * x3 + px[2] * x0 + px[3] * x1;
                ty = py[0] * x4 + py[1] * x3 + py[2] * x0 + py[3] * x1;
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
                //	- Result: 4*ADD + 2*ASSIGNMENT + 26*MUL + 2*POW
                //	- Subexpressions: 60*ADD + DIV + 191*MUL + 11*POW + 56*SUB
                T x0 = px[2] * py[1];
                T x1 = px[3] * py[0];
                T x2 = px[0] * py[3];
                T x3 = px[1] * py[2];
                T x4 = px[0] * py[2];
                T x5 = px[1] * py[3];
                T x6 = px[2] * py[0];
                T x7 = px[3] * py[1];
                T x8 = x4 + x5 - x6 - x7;
                T x9 = 2 * x4;
                T x10 = px[1] * ty;
                T x11 = 2 * x10;
                T x12 = px[0] * ty;
                T x13 = 2 * x12;
                T x14 = 4 * x2;
                T x15 = px[2] * py[3];
                T x16 = px[2] * ty;
                T x17 = px[3] * ty;
                T x18 = py[0] * tx;
                T x19 = py[1] * tx;
                T x20 = py[3] * tx;
                T x21 = 2 * tx;
                T x22 = ty * x21;
                T x23 = 2 * x2;
                T x24 = tx * ty;
                T x25 = 2 * x5;
                T x26 = 4 * x3;
                T x27 = 4 * x1;
                T x28 = px[3] * py[2];
                T x29 = 2 * x17;
                T x30 = py[2] * tx;
                T x31 = 2 * x30;
                T x32 = 2 * x1;
                T x33 = 6 * x1;
                T x34 = 2 * px[3];
                T x35 = py[3] * x34;
                T x36 = 2 * x28;
                T x37 = 2 * x20;
                T x38 = 2 * px[2];
                T x39 = pow(px[0], 2);
                T x40 = pow(py[2], 2);
                T x41 = pow(ty, 2);
                T x42 = pow(px[1], 2);
                T x43 = pow(py[3], 2);
                T x44 = pow(px[2], 2);
                T x45 = pow(py[0], 2);
                T x46 = pow(px[3], 2);
                T x47 = pow(py[1], 2);
                T x48 = pow(tx, 2);
                T x49 = 2 * px[1];
                T x50 = px[0] * x41;
                T x51 = 4 * x43;
                T x52 = 2 * x40;
                T x53 = px[3] * x52;
                T x54 = px[0] * tx;
                T x55 = 2 * ty;
                T x56 = py[2] * x55;
                T x57 = px[2] * x49;
                T x58 = px[2] * x21;
                T x59 = 2 * x44;
                T x60 = py[3] * x59;
                T x61 = py[0] * ty;
                T x62 = 4 * x46;
                T x63 = 2 * py[1];
                T x64 = py[2] * x63;
                T x65 = py[0] * x48;
                T x66 = 2 * py[3];
                T x67 =
                    sqrt(px[0] * px[2] * x51 + px[0] * x53 + px[1] * x21 * x43 - px[2] * py[2] * x32 -
                         px[2] * x34 * x41 - 4 * px[3] * py[3] * x4 - px[3] * x21 * x47 + px[3] * x41 * x49 +
                         py[0] * py[2] * x62 + py[0] * x60 + py[1] * x46 * x55 + py[2] * x20 * x38 - py[2] * x48 * x66 +
                         2 * py[2] * x65 - py[3] * x42 * x55 + py[3] * x48 * x63 - tx * x53 - ty * x60 - 4 * x0 * x12 +
                         x0 * x14 + 2 * x0 * x18 + x0 * x22 - x0 * x29 - x0 * x32 + x0 * x35 - x0 * x37 + x1 * x26 -
                         x10 * x27 + x11 * x18 - x11 * x19 + x11 * x4 + x11 * x6 + x11 * x7 - x13 * x18 + x13 * x19 +
                         x13 * x5 + x13 * x6 + x13 * x7 + x14 * x17 - x14 * x19 - x15 * x27 - x15 * x36 - x15 * x9 -
                         6 * x16 * x2 - x16 * x31 + x16 * x33 + x16 * x35 + x16 * x36 + x16 * x37 + 4 * x16 * x5 +
                         x16 * x9 - x17 * x37 + x18 * x25 - x18 * x26 + x18 * x9 + x19 * x25 + x19 * x32 + x19 * x9 +
                         x20 * x27 + 6 * x20 * x4 + x22 * x28 + x22 * x3 - x22 * x4 - x22 * x5 - x22 * x6 - x22 * x7 +
                         x23 * x24 + x24 * x32 + x25 * x28 - x25 * x6 - x25 * x7 + x28 * x37 - x29 * x3 - x29 * x5 -
                         x30 * x33 + 4 * x30 * x7 - x31 * x5 + x31 * x6 - x34 * x50 - x37 * x7 + x38 * x50 + x39 * x40 +
                         x39 * x41 - x39 * x56 + x40 * x46 + x40 * x48 + x41 * x42 + x41 * x44 + x41 * x46 - x41 * x57 +
                         x42 * x43 + x43 * x44 + x43 * x48 - x43 * x57 - x43 * x58 + x44 * x45 + x45 * x48 - x45 * x58 +
                         x46 * x47 - x46 * x56 - x46 * x64 + x47 * x48 - x48 * x64 - x49 * x50 - x5 * x9 - x51 * x54 -
                         x52 * x54 - x59 * x61 - x6 * x9 - x61 * x62 - x63 * x65 - x65 * x66 - x7 * x9);
                T x68 = -x4 + x6;
                T x69 = px[1] * py[0];
                x = (1.0 / 2.0) *
                    (-x10 + x12 - x15 + x16 - x17 - x18 + x19 + x20 - x23 + x28 - x30 + x32 + x5 - x67 - x68 - x7) /
                    (x0 + x1 - x2 - x3 + x8);
                y = (1.0 / 2.0) *
                    (2 * px[0] * py[1] + x10 - x12 + x15 - x16 + x17 + x18 - x19 - x20 + x23 - x28 + x30 - x32 - x67 -
                     2 * x69 - x8) /
                    (px[0] * py[1] - x1 + x2 + x68 - x69);
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
                //	- Result: 4*ADD + 8*ASSIGNMENT + 12*MUL
                //	- Subexpressions: 15*ADD + DIV + 38*MUL + 2*NEG + 22*SUB
                T x0 = x - 1;
                T x1 = -px[0] * x + px[0] + px[1] * x - px[2] * x + px[3] * x;
                T x2 = y - 1;
                T x3 = py[1] * x;
                T x4 = py[3] * x;
                T x5 = py[0] * x;
                T x6 = py[2] * x;
                T x7 = py[0] + x3 + x4 - x5 - x6;
                T x8 = -x2 * x7;
                T x9 = px[1] * py[0];
                T x10 = px[0] * py[1];
                T x11 = py[3] * y;
                T x12 = py[0] * y;
                T x13 = py[2] * y;
                T x14 = 1.0 /
                        (px[0] * x11 - px[0] * x13 - px[0] * x4 + px[0] * x6 + px[1] * x4 - px[1] * x6 + px[2] * x12 +
                         px[2] * x3 - px[2] * x5 - px[3] * x12 - px[3] * x3 + px[3] * x5 + x10 * y - x10 - x9 * y + x9);
                T x15 = px[0] * y - px[0] - px[1] * y + px[1] + px[2] * y - px[3] * y;
                T x16 = -py[0] - py[1] * y + py[1] - x11 + x12 + x13;
                T x17 = -x16 * x2;
                T x18 = x * x1;
                T x19 = x * x15;
                T x20 = x18 - x7 * y;
                T x21 = -x16 * y + x19;
                gx[0] = x14 * (-x0 * x1 - x8);
                gy[0] = x14 * (-x0 * x15 - x17);
                gx[1] = x14 * (x18 + x8);
                gy[1] = x14 * (x17 + x19);
                gx[2] = -x14 * x20;
                gy[2] = -x14 * x21;
                gx[3] = x14 * x20;
                gy[3] = x14 * x21;
            }

            UTOPIA_FUNCTION static void value(const T x, const T y, Result *UTOPIA_RESTRICT f) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 2*ADD + 4*ASSIGNMENT + 4*MUL
                //	- Subexpressions: MUL + SUB
                T x0 = 1 - y;
                T x1 = x * y;
                f[0] = x0 * (1 - x);
                f[1] = x * x0;
                f[2] = x1;
                f[3] = 1 - x1;
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
                //	- Result: 5*ADD + 13*ASSIGNMENT + 16*MUL
                //	- Subexpressions: 15*ADD + DIV + 47*MUL + 3*NEG + 22*SUB
                T x0 = x - 1;
                T x1 = y - 1;
                T x2 = -x1;
                T x3 = x * y;
                T x4 = py[2] * x;
                T x5 = py[3] * y;
                T x6 = py[3] * x;
                T x7 = py[0] * y;
                T x8 = py[1] * x;
                T x9 = py[0] * x;
                T x10 = px[0] * py[1] * y - px[0] * py[1] - px[0] * py[2] * y - px[0] * py[3] * x + px[0] * x4 +
                        px[0] * x5 - px[1] * py[0] * y + px[1] * py[0] - px[1] * py[2] * x + px[1] * x6 -
                        px[2] * py[0] * x + px[2] * x7 + px[2] * x8 - px[3] * py[0] * y - px[3] * py[1] * x +
                        px[3] * x9;
                T x11 = -px[0] * x + px[0] + px[1] * x - px[2] * x + px[3] * x;
                T x12 = py[0] - x4 + x6 + x8 - x9;
                T x13 = -x1 * x12;
                T x14 = 1.0 / x10;
                T x15 = px[0] * y - px[0] - px[1] * y + px[1] + px[2] * y - px[3] * y;
                T x16 = -py[0] - py[1] * y + py[1] + py[2] * y - x5 + x7;
                T x17 = -x1 * x16;
                T x18 = x * x11;
                T x19 = x * x15;
                T x20 = -x12 * y + x18;
                T x21 = -x16 * y + x19;
                f[0] = -x0 * x2;
                f[1] = x * x2;
                f[2] = x3;
                f[3] = 1 - x3;
                measure_value = -x10;
                gx[0] = x14 * (-x0 * x11 - x13);
                gy[0] = x14 * (-x0 * x15 - x17);
                gx[1] = x14 * (x13 + x18);
                gy[1] = x14 * (x17 + x19);
                gx[2] = -x14 * x20;
                gy[2] = -x14 * x21;
                gx[3] = x14 * x20;
                gy[3] = x14 * x21;
            }
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_FE_Quad4_2_IMPL_hpp
