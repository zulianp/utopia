#ifndef UTOPIA_TPL_FE_AxisAlignedHex8_3_IMPL_hpp
#define UTOPIA_TPL_FE_AxisAlignedHex8_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

// #include "utopia_fe_AxisAlignedHex8.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of AxisAlignedHex8 for dimension 3
         */
        template <typename T, typename GeoT>
        class AxisAlignedHex8 {
        public:
            static constexpr int Dim = 3;
            static constexpr int NNodes = 8;
            static constexpr int Order = 1;

            using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "AxisAlignedHex8"; }

            UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

            UTOPIA_INLINE_FUNCTION static constexpr int n_nodes() { return NNodes; }

            UTOPIA_INLINE_FUNCTION static constexpr int order() { return Order; }

            UTOPIA_FUNCTION static constexpr Result measure(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Input quadrature point
                const T x,
                const T y,
                const T z) {
                T measure_value;

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                // FLOATING POINT OPS!
                //	- Result: 3*ADD + ASSIGNMENT + 4*MUL
                //	- Subexpressions: 0
                measure_value = (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]);
                return measure_value;
            }

            UTOPIA_FUNCTION static void jacobian(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Input quadrature point
                const T x,
                const T y,
                const T z,
                GeoT *UTOPIA_RESTRICT J) {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                // FLOATING POINT OPS!
                //	- Result: 3*ADD + 9*ASSIGNMENT + 3*MUL
                //	- Subexpressions: 0
                J[0] = -px[0] + px[6];
                J[1] = 0;
                J[2] = 0;
                J[3] = 0;
                J[4] = -py[0] + py[6];
                J[5] = 0;
                J[6] = 0;
                J[7] = 0;
                J[8] = -pz[0] + pz[6];
            }

            UTOPIA_FUNCTION static void jacobian_inverse(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Input quadrature point
                const T x,
                const T y,
                const T z,
                GeoT *UTOPIA_RESTRICT J_inv) {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                // FLOATING POINT OPS!
                //	- Result: 3*ADD + 9*ASSIGNMENT + 6*MUL + 3*POW
                //	- Subexpressions: 0
                J_inv[0] = -1 / (px[0] - px[6]);
                J_inv[1] = 0;
                J_inv[2] = 0;
                J_inv[3] = 0;
                J_inv[4] = -1 / (py[0] - py[6]);
                J_inv[5] = 0;
                J_inv[6] = 0;
                J_inv[7] = 0;
                J_inv[8] = -1 / (pz[0] - pz[6]);
            }

            UTOPIA_FUNCTION static void transform(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Input quadrature point
                const T x,
                const T y,
                const T z,
                GeoT &tx,
                GeoT &ty,
                GeoT &tz) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 6*ADD + 3*ASSIGNMENT + 6*MUL
                //	- Subexpressions: 0
                tx = px[0] - x * (px[0] - px[6]);
                ty = py[0] - y * (py[0] - py[6]);
                tz = pz[0] - z * (pz[0] - pz[6]);
            }

            UTOPIA_FUNCTION static void inverse_transform(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Input quadrature point
                const T tx,
                const T ty,
                const T tz,
                GeoT &x,
                GeoT &y,
                GeoT &z) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 6*ADD + 3*ASSIGNMENT + 9*MUL + 3*POW
                //	- Subexpressions: 0
                x = (px[0] - tx) / (px[0] - px[6]);
                y = (py[0] - ty) / (py[0] - py[6]);
                z = (pz[0] - tz) / (pz[0] - pz[6]);
            }

            UTOPIA_FUNCTION static void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Input quadrature point
                const T x,
                const T y,
                const T z,
                // Output
                Result *UTOPIA_RESTRICT gx,
                Result *UTOPIA_RESTRICT gy,
                Result *UTOPIA_RESTRICT gz) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 24*ASSIGNMENT + 24*MUL
                //	- Subexpressions: 3*ADD + 8*DIV + 20*MUL + 9*SUB
                T x0 = y - 1.0;
                T x1 = 1.0 / (px[0] - px[6]);
                T x2 = z - 1.0;
                T x3 = (1.0 / 8.0) * x2;
                T x4 = x1 * x3;
                T x5 = x - 1.0;
                T x6 = 1.0 / (py[0] - py[6]);
                T x7 = x3 * x6;
                T x8 = (1.0 / 8.0) / (pz[0] - pz[6]);
                T x9 = x0 * x8;
                T x10 = 1.0 / (8 * px[0] - 8 * px[6]);
                T x11 = x10 * x2;
                T x12 = x + 1.0;
                T x13 = 1.0 / (8 * py[0] - 8 * py[6]);
                T x14 = x13 * x2;
                T x15 = 1.0 / (8 * pz[0] - 8 * pz[6]);
                T x16 = x0 * x15;
                T x17 = y + 1.0;
                T x18 = x17 * x8;
                T x19 = x15 * x17;
                T x20 = z + 1.0;
                T x21 = x10 * x20;
                T x22 = x13 * x20;
                T x23 = (1.0 / 8.0) * x20;
                T x24 = x1 * x23;
                T x25 = x23 * x6;
                gx[0] = x0 * x4;
                gy[0] = x5 * x7;
                gz[0] = x5 * x9;
                gx[1] = -x0 * x11;
                gy[1] = -x12 * x14;
                gz[1] = -x12 * x16;
                gx[2] = x17 * x4;
                gy[2] = x12 * x7;
                gz[2] = x12 * x18;
                gx[3] = -x11 * x17;
                gy[3] = -x14 * x5;
                gz[3] = -x19 * x5;
                gx[4] = -x0 * x21;
                gy[4] = -x22 * x5;
                gz[4] = -x16 * x5;
                gx[5] = x0 * x24;
                gy[5] = x12 * x25;
                gz[5] = x12 * x9;
                gx[6] = -x17 * x21;
                gy[6] = -x12 * x22;
                gz[6] = -x12 * x19;
                gx[7] = x17 * x24;
                gy[7] = x25 * x5;
                gz[7] = x18 * x5;
            }

            UTOPIA_FUNCTION static void value(const T x, const T y, const T z, Result *UTOPIA_RESTRICT f) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 8*ASSIGNMENT + 8*MUL
                //	- Subexpressions: 3*ADD + 3*DIV + 4*MUL + 3*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * z;
                T x5 = 0.5 - x4;
                T x6 = x3 * x5;
                T x7 = x0 + 0.5;
                T x8 = x2 + 0.5;
                T x9 = x5 * x8;
                T x10 = x4 + 0.5;
                T x11 = x10 * x3;
                T x12 = x10 * x8;
                f[0] = x1 * x6;
                f[1] = x6 * x7;
                f[2] = x7 * x9;
                f[3] = x1 * x9;
                f[4] = x1 * x11;
                f[5] = x11 * x7;
                f[6] = x12 * x7;
                f[7] = x1 * x12;
            }

            UTOPIA_FUNCTION static void eval(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Input quadrature point
                const T x,
                const T y,
                const T z,
                // Output
                Result *UTOPIA_RESTRICT f,
                Result *UTOPIA_RESTRICT gx,
                Result *UTOPIA_RESTRICT gy,
                Result *UTOPIA_RESTRICT gz,
                T &measure_value) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 33*ASSIGNMENT + 33*MUL
                //	- Subexpressions: 6*ADD + 11*DIV + 23*MUL + 12*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * z;
                T x5 = 0.5 - x4;
                T x6 = x3 * x5;
                T x7 = x0 + 0.5;
                T x8 = x2 + 0.5;
                T x9 = x5 * x8;
                T x10 = x4 + 0.5;
                T x11 = x10 * x3;
                T x12 = x10 * x8;
                T x13 = px[0] - px[6];
                T x14 = py[0] - py[6];
                T x15 = pz[0] - pz[6];
                T x16 = y - 1.0;
                T x17 = 1.0 / x13;
                T x18 = z - 1.0;
                T x19 = (1.0 / 8.0) * x18;
                T x20 = x17 * x19;
                T x21 = x - 1.0;
                T x22 = 1.0 / x14;
                T x23 = x19 * x22;
                T x24 = (1.0 / 8.0) / x15;
                T x25 = x16 * x24;
                T x26 = 1.0 / (8 * px[0] - 8 * px[6]);
                T x27 = x18 * x26;
                T x28 = x + 1.0;
                T x29 = 1.0 / (8 * py[0] - 8 * py[6]);
                T x30 = x18 * x29;
                T x31 = 1.0 / (8 * pz[0] - 8 * pz[6]);
                T x32 = x16 * x31;
                T x33 = y + 1.0;
                T x34 = x24 * x33;
                T x35 = x31 * x33;
                T x36 = z + 1.0;
                T x37 = x26 * x36;
                T x38 = x29 * x36;
                T x39 = (1.0 / 8.0) * x36;
                T x40 = x17 * x39;
                T x41 = x22 * x39;
                f[0] = x1 * x6;
                f[1] = x6 * x7;
                f[2] = x7 * x9;
                f[3] = x1 * x9;
                f[4] = x1 * x11;
                f[5] = x11 * x7;
                f[6] = x12 * x7;
                f[7] = x1 * x12;
                measure_value = -x13 * x14 * x15;
                gx[0] = x16 * x20;
                gy[0] = x21 * x23;
                gz[0] = x21 * x25;
                gx[1] = -x16 * x27;
                gy[1] = -x28 * x30;
                gz[1] = -x28 * x32;
                gx[2] = x20 * x33;
                gy[2] = x23 * x28;
                gz[2] = x28 * x34;
                gx[3] = -x27 * x33;
                gy[3] = -x21 * x30;
                gz[3] = -x21 * x35;
                gx[4] = -x16 * x37;
                gy[4] = -x21 * x38;
                gz[4] = -x21 * x32;
                gx[5] = x16 * x40;
                gy[5] = x28 * x41;
                gz[5] = x25 * x28;
                gx[6] = -x33 * x37;
                gy[6] = -x28 * x38;
                gz[6] = -x28 * x35;
                gx[7] = x33 * x40;
                gy[7] = x21 * x41;
                gz[7] = x21 * x34;
            }
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_FE_AxisAlignedHex8_3_IMPL_hpp
