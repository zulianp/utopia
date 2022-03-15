#ifndef UTOPIA_TPL_FE_Tet4_3_IMPL_hpp
#define UTOPIA_TPL_FE_Tet4_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

// #include "utopia_fe_Tet4.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Tet4 for dimension 3
         */
        template <typename T, typename GeoT = T>
        class Tet4 {
        public:
            static constexpr int Dim = 3;
            static constexpr int NNodes = 4;
            static constexpr int Order = 1;

            using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Tet4"; }

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
                // FLOATING POINT OPS!
                //	- Result: ADD + ASSIGNMENT + 6*MUL
                //	- Subexpressions: 9*SUB
                T x0 = -px[0] + px[1];
                T x1 = -py[0] + py[2];
                T x2 = -pz[0] + pz[3];
                T x3 = -px[0] + px[2];
                T x4 = -py[0] + py[3];
                T x5 = -pz[0] + pz[1];
                T x6 = -px[0] + px[3];
                T x7 = -py[0] + py[1];
                T x8 = -pz[0] + pz[2];
                measure_value = x0 * x1 * x2 - x0 * x4 * x8 - x1 * x5 * x6 - x2 * x3 * x7 + x3 * x4 * x5 + x6 * x7 * x8;
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
                // FLOATING POINT OPS!
                //	- Result: 9*ADD + 9*ASSIGNMENT + 9*MUL
                //	- Subexpressions: 0
                J[0] = -px[0] + px[1];
                J[1] = -py[0] + py[1];
                J[2] = -pz[0] + pz[1];
                J[3] = -px[0] + px[2];
                J[4] = -py[0] + py[2];
                J[5] = -pz[0] + pz[2];
                J[6] = -px[0] + px[3];
                J[7] = -py[0] + py[3];
                J[8] = -pz[0] + pz[3];
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
                // FLOATING POINT OPS!
                //	- Result: 9*ADD + 9*ASSIGNMENT + 25*MUL
                //	- Subexpressions: 2*ADD + DIV + 12*MUL + 12*SUB
                T x0 = -py[0] + py[2];
                T x1 = -pz[0] + pz[3];
                T x2 = x0 * x1;
                T x3 = -py[0] + py[3];
                T x4 = -pz[0] + pz[2];
                T x5 = x3 * x4;
                T x6 = -px[0] + px[1];
                T x7 = -px[0] + px[2];
                T x8 = -pz[0] + pz[1];
                T x9 = -px[0] + px[3];
                T x10 = -py[0] + py[1];
                T x11 = x10 * x4;
                T x12 = x1 * x10;
                T x13 = x0 * x8;
                T x14 = 1.0 / (x11 * x9 - x12 * x7 - x13 * x9 + x2 * x6 + x3 * x7 * x8 - x5 * x6);
                J_inv[0] = x14 * (x2 - x5);
                J_inv[1] = x14 * (-x12 + x3 * x8);
                J_inv[2] = x14 * (x11 - x13);
                J_inv[3] = x14 * (-x1 * x7 + x4 * x9);
                J_inv[4] = x14 * (x1 * x6 - x8 * x9);
                J_inv[5] = x14 * (-x4 * x6 + x7 * x8);
                J_inv[6] = x14 * (-x0 * x9 + x3 * x7);
                J_inv[7] = x14 * (x10 * x9 - x3 * x6);
                J_inv[8] = x14 * (x0 * x6 - x10 * x7);
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
                //	- Result: 3*ADD + 3*ASSIGNMENT + 12*MUL
                //	- Subexpressions: 3*SUB
                T x0 = -x - y - z + 1;
                tx = px[0] * x0 + px[1] * x + px[2] * y + px[3] * z;
                ty = py[0] * x0 + py[1] * x + py[2] * y + py[3] * z;
                tz = pz[0] * x0 + pz[1] * x + pz[2] * y + pz[3] * z;
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
                //	- Result: 3*ADD + 3*ASSIGNMENT + 35*MUL
                //	- Subexpressions: 15*ADD + DIV + 62*MUL + 18*SUB
                T x0 = py[2] * tx;
                T x1 = pz[0] * x0;
                T x2 = py[3] * tx;
                T x3 = px[0] * ty;
                T x4 = pz[2] * x3;
                T x5 = px[2] * ty;
                T x6 = px[2] * py[0];
                T x7 = tz * x6;
                T x8 = px[3] * py[2];
                T x9 = py[0] * tx;
                T x10 = pz[2] * x9;
                T x11 = pz[0] * x5;
                T x12 = px[3] * ty;
                T x13 = px[0] * py[2];
                T x14 = tz * x13;
                T x15 = px[2] * py[3];
                T x16 = px[3] * py[0];
                T x17 = px[0] * py[3];
                T x18 = pz[0] * x15 - pz[0] * x8 + pz[2] * x16 - pz[2] * x17 + pz[3] * x13 - pz[3] * x6;
                T x19 = pz[0] * x12 - pz[0] * x2 - pz[3] * x3 + pz[3] * x9 - tz * x16 + tz * x17;
                T x20 = px[1] * py[0];
                T x21 = px[1] * py[3];
                T x22 = px[2] * py[1];
                T x23 = px[3] * py[1];
                T x24 = px[0] * py[1];
                T x25 = pz[3] * x24;
                T x26 = px[1] * py[2];
                T x27 = pz[0] * x21;
                T x28 = pz[1] * x16;
                T x29 = -pz[0] * x22 + pz[0] * x26 - pz[1] * x13 + pz[1] * x6 - pz[2] * x20 + pz[2] * x24;
                T x30 = 1.0 / (pz[0] * x23 - pz[1] * x15 + pz[1] * x17 + pz[1] * x8 + pz[2] * x21 - pz[2] * x23 +
                               pz[3] * x20 + pz[3] * x22 - pz[3] * x26 + x18 - x25 - x27 - x28 + x29);
                T x31 = px[1] * ty;
                T x32 = py[1] * tx;
                T x33 = -pz[0] * x31 + pz[0] * x32 + pz[1] * x3 - pz[1] * x9 + tz * x20 - tz * x24;
                x = x30 * (-pz[2] * x12 + pz[2] * x2 - pz[3] * x0 + pz[3] * x5 - tz * x15 + tz * x8 + x1 - x10 - x11 -
                           x14 + x18 + x19 + x4 + x7);
                y = x30 * (px[0] * py[3] * pz[1] + px[1] * py[0] * pz[3] + px[1] * py[3] * tz + px[3] * py[1] * pz[0] +
                           px[3] * pz[1] * ty + py[1] * pz[3] * tx - pz[1] * x2 - pz[3] * x31 - tz * x23 - x19 - x25 -
                           x27 - x28 - x33);
                z = x30 * (pz[1] * x0 - pz[1] * x5 + pz[2] * x31 - pz[2] * x32 + tz * x22 - tz * x26 - x1 + x10 + x11 +
                           x14 + x29 + x33 - x4 - x7);
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
                //	- Result: 3*ADD + 12*ASSIGNMENT + 9*MUL
                //	- Subexpressions: 2*ADD + DIV + 34*MUL + 21*SUB
                T x0 = -py[0] + py[1];
                T x1 = -pz[0] + pz[2];
                T x2 = x0 * x1;
                T x3 = -py[0] + py[2];
                T x4 = -pz[0] + pz[1];
                T x5 = x3 * x4;
                T x6 = -px[0] + px[1];
                T x7 = -pz[0] + pz[3];
                T x8 = x3 * x7;
                T x9 = -px[0] + px[2];
                T x10 = -py[0] + py[3];
                T x11 = -px[0] + px[3];
                T x12 = x1 * x10;
                T x13 = x0 * x7;
                T x14 = 1.0 / (x10 * x4 * x9 + x11 * x2 - x11 * x5 - x12 * x6 - x13 * x9 + x6 * x8);
                T x15 = x14 * (x2 - x5);
                T x16 = x14 * (x10 * x4 - x13);
                T x17 = x14 * (-x12 + x8);
                T x18 = x14 * (-x1 * x6 + x4 * x9);
                T x19 = x14 * (-x11 * x4 + x6 * x7);
                T x20 = x14 * (x1 * x11 - x7 * x9);
                T x21 = x14 * (-x0 * x9 + x3 * x6);
                T x22 = x14 * (x0 * x11 - x10 * x6);
                T x23 = x14 * (x10 * x9 - x11 * x3);
                gx[0] = -x15 - x16 - x17;
                gy[0] = -x18 - x19 - x20;
                gz[0] = -x21 - x22 - x23;
                gx[1] = x17;
                gy[1] = x20;
                gz[1] = x23;
                gx[2] = x16;
                gy[2] = x19;
                gz[2] = x22;
                gx[3] = x15;
                gy[3] = x18;
                gz[3] = x21;
            }

            UTOPIA_FUNCTION static void value(const T x, const T y, const T z, Result *UTOPIA_RESTRICT f) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: ADD + 4*ASSIGNMENT + 3*MUL
                //	- Subexpressions: 0
                f[0] = -x - y - z + 1;
                f[1] = x;
                f[2] = y;
                f[3] = z;
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
                //	- Result: 4*ADD + 17*ASSIGNMENT + 12*MUL
                //	- Subexpressions: 2*ADD + DIV + 34*MUL + 21*SUB
                T x0 = -px[0] + px[1];
                T x1 = -py[0] + py[2];
                T x2 = -pz[0] + pz[3];
                T x3 = x1 * x2;
                T x4 = -px[0] + px[2];
                T x5 = -py[0] + py[3];
                T x6 = -pz[0] + pz[1];
                T x7 = -px[0] + px[3];
                T x8 = -py[0] + py[1];
                T x9 = -pz[0] + pz[2];
                T x10 = x8 * x9;
                T x11 = x5 * x9;
                T x12 = x2 * x8;
                T x13 = x1 * x6;
                T x14 = -x0 * x11 + x0 * x3 + x10 * x7 - x12 * x4 - x13 * x7 + x4 * x5 * x6;
                T x15 = 1.0 / x14;
                T x16 = x15 * (x10 - x13);
                T x17 = x15 * (-x12 + x5 * x6);
                T x18 = x15 * (-x11 + x3);
                T x19 = x15 * (-x0 * x9 + x4 * x6);
                T x20 = x15 * (x0 * x2 - x6 * x7);
                T x21 = x15 * (-x2 * x4 + x7 * x9);
                T x22 = x15 * (x0 * x1 - x4 * x8);
                T x23 = x15 * (-x0 * x5 + x7 * x8);
                T x24 = x15 * (-x1 * x7 + x4 * x5);
                f[0] = -x - y - z + 1;
                f[1] = x;
                f[2] = y;
                f[3] = z;
                measure_value = x14;
                gx[0] = -x16 - x17 - x18;
                gy[0] = -x19 - x20 - x21;
                gz[0] = -x22 - x23 - x24;
                gx[1] = x18;
                gy[1] = x21;
                gz[1] = x24;
                gx[2] = x17;
                gy[2] = x20;
                gz[2] = x23;
                gx[3] = x16;
                gy[3] = x19;
                gz[3] = x22;
            }
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_FE_Tet4_3_IMPL_hpp
