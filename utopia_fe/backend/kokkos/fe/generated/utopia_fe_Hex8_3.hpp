#ifndef UTOPIA_TPL_FE_Hex8_3_IMPL_hpp
#define UTOPIA_TPL_FE_Hex8_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

// #include "utopia_fe_Hex8.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Hex8 for dimension 3
         */
        template <typename T, typename GeoT>
        class Hex8 {
        public:
            static constexpr int Dim = 3;
            static constexpr int NNodes = 8;
            static constexpr int Order = 1;

            using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Hex8"; }

            UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

            UTOPIA_INLINE_FUNCTION static constexpr int n_nodes() { return NNodes; }

            UTOPIA_INLINE_FUNCTION static constexpr int order() { return Order; }

            UTOPIA_FUNCTION static Result measure(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Input quadrature point
                const T x,
                const T y,
                const T z) {
                T measure_value = 0;
                // FLOATING POINT OPS!
                //	- Result: ADD + ASSIGNMENT + 6*MUL
                //	- Subexpressions: 30*ADD + 43*DIV + 120*MUL + 39*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * x3;
                T x5 = x1 * x4;
                T x6 = x0 + 0.5;
                T x7 = x4 * x6;
                T x8 = x2 + 0.5;
                T x9 = (1.0 / 2.0) * x8;
                T x10 = x6 * x9;
                T x11 = x1 * x9;
                T x12 = -px[0] * x5 - px[1] * x7 - px[2] * x10 - px[3] * x11 + (1.0 / 2.0) * px[4] * x1 * x3 +
                        (1.0 / 2.0) * px[5] * x3 * x6 + (1.0 / 2.0) * px[6] * x6 * x8 + (1.0 / 2.0) * px[7] * x1 * x8;
                T x13 = (1.0 / 2.0) * z;
                T x14 = 0.5 - x13;
                T x15 = (1.0 / 2.0) * x14;
                T x16 = x1 * x15;
                T x17 = x15 * x6;
                T x18 = x13 + 0.5;
                T x19 = (1.0 / 2.0) * x18;
                T x20 = x1 * x19;
                T x21 = x19 * x6;
                T x22 = -py[0] * x16 - py[1] * x17 + (1.0 / 2.0) * py[2] * x14 * x6 + (1.0 / 2.0) * py[3] * x1 * x14 -
                        py[4] * x20 - py[5] * x21 + (1.0 / 2.0) * py[6] * x18 * x6 + (1.0 / 2.0) * py[7] * x1 * x18;
                T x23 = x15 * x3;
                T x24 = x15 * x8;
                T x25 = x19 * x3;
                T x26 = x19 * x8;
                T x27 = -pz[0] * x23 + (1.0 / 2.0) * pz[1] * x14 * x3 + (1.0 / 2.0) * pz[2] * x14 * x8 - pz[3] * x24 -
                        pz[4] * x25 + (1.0 / 2.0) * pz[5] * x18 * x3 + (1.0 / 2.0) * pz[6] * x18 * x8 - pz[7] * x26;
                T x28 = -px[0] * x16 - px[1] * x17 + (1.0 / 2.0) * px[2] * x14 * x6 + (1.0 / 2.0) * px[3] * x1 * x14 -
                        px[4] * x20 - px[5] * x21 + (1.0 / 2.0) * px[6] * x18 * x6 + (1.0 / 2.0) * px[7] * x1 * x18;
                T x29 = -py[0] * x23 + (1.0 / 2.0) * py[1] * x14 * x3 + (1.0 / 2.0) * py[2] * x14 * x8 - py[3] * x24 -
                        py[4] * x25 + (1.0 / 2.0) * py[5] * x18 * x3 + (1.0 / 2.0) * py[6] * x18 * x8 - py[7] * x26;
                T x30 = -pz[0] * x5 - pz[1] * x7 - pz[2] * x10 - pz[3] * x11 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                        (1.0 / 2.0) * pz[5] * x3 * x6 + (1.0 / 2.0) * pz[6] * x6 * x8 + (1.0 / 2.0) * pz[7] * x1 * x8;
                T x31 = -px[0] * x23 + (1.0 / 2.0) * px[1] * x14 * x3 + (1.0 / 2.0) * px[2] * x14 * x8 - px[3] * x24 -
                        px[4] * x25 + (1.0 / 2.0) * px[5] * x18 * x3 + (1.0 / 2.0) * px[6] * x18 * x8 - px[7] * x26;
                T x32 = -py[0] * x5 - py[1] * x7 - py[2] * x10 - py[3] * x11 + (1.0 / 2.0) * py[4] * x1 * x3 +
                        (1.0 / 2.0) * py[5] * x3 * x6 + (1.0 / 2.0) * py[6] * x6 * x8 + (1.0 / 2.0) * py[7] * x1 * x8;
                T x33 = -pz[0] * x16 - pz[1] * x17 + (1.0 / 2.0) * pz[2] * x14 * x6 + (1.0 / 2.0) * pz[3] * x1 * x14 -
                        pz[4] * x20 - pz[5] * x21 + (1.0 / 2.0) * pz[6] * x18 * x6 + (1.0 / 2.0) * pz[7] * x1 * x18;
                measure_value = -x12 * x22 * x27 + x12 * x29 * x33 + x22 * x30 * x31 + x27 * x28 * x32 -
                                x28 * x29 * x30 - x31 * x32 * x33;
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
                //	- Result: 9*ADD + 9*ASSIGNMENT + 72*MUL
                //	- Subexpressions: 3*ADD + 7*DIV + 12*MUL + 3*SUB
                T x0 = (1.0 / 2.0) * y;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * z;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * x3;
                T x5 = x1 * x4;
                T x6 = x0 + 0.5;
                T x7 = x4 * x6;
                T x8 = x2 + 0.5;
                T x9 = (1.0 / 2.0) * x8;
                T x10 = x1 * x9;
                T x11 = x6 * x9;
                T x12 = (1.0 / 2.0) * x;
                T x13 = 0.5 - x12;
                T x14 = x13 * x4;
                T x15 = x12 + 0.5;
                T x16 = x15 * x4;
                T x17 = x13 * x9;
                T x18 = x15 * x9;
                T x19 = (1.0 / 2.0) * x1;
                T x20 = x13 * x19;
                T x21 = x15 * x19;
                T x22 = (1.0 / 2.0) * x6;
                T x23 = x15 * x22;
                T x24 = x13 * x22;
                J[0] = -px[0] * x5 + (1.0 / 2.0) * px[1] * x1 * x3 + (1.0 / 2.0) * px[2] * x3 * x6 - px[3] * x7 -
                       px[4] * x10 + (1.0 / 2.0) * px[5] * x1 * x8 + (1.0 / 2.0) * px[6] * x6 * x8 - px[7] * x11;
                J[1] = -px[0] * x14 - px[1] * x16 + (1.0 / 2.0) * px[2] * x15 * x3 + (1.0 / 2.0) * px[3] * x13 * x3 -
                       px[4] * x17 - px[5] * x18 + (1.0 / 2.0) * px[6] * x15 * x8 + (1.0 / 2.0) * px[7] * x13 * x8;
                J[2] = -px[0] * x20 - px[1] * x21 - px[2] * x23 - px[3] * x24 + (1.0 / 2.0) * px[4] * x1 * x13 +
                       (1.0 / 2.0) * px[5] * x1 * x15 + (1.0 / 2.0) * px[6] * x15 * x6 + (1.0 / 2.0) * px[7] * x13 * x6;
                J[3] = -py[0] * x5 + (1.0 / 2.0) * py[1] * x1 * x3 + (1.0 / 2.0) * py[2] * x3 * x6 - py[3] * x7 -
                       py[4] * x10 + (1.0 / 2.0) * py[5] * x1 * x8 + (1.0 / 2.0) * py[6] * x6 * x8 - py[7] * x11;
                J[4] = -py[0] * x14 - py[1] * x16 + (1.0 / 2.0) * py[2] * x15 * x3 + (1.0 / 2.0) * py[3] * x13 * x3 -
                       py[4] * x17 - py[5] * x18 + (1.0 / 2.0) * py[6] * x15 * x8 + (1.0 / 2.0) * py[7] * x13 * x8;
                J[5] = -py[0] * x20 - py[1] * x21 - py[2] * x23 - py[3] * x24 + (1.0 / 2.0) * py[4] * x1 * x13 +
                       (1.0 / 2.0) * py[5] * x1 * x15 + (1.0 / 2.0) * py[6] * x15 * x6 + (1.0 / 2.0) * py[7] * x13 * x6;
                J[6] = -pz[0] * x5 + (1.0 / 2.0) * pz[1] * x1 * x3 + (1.0 / 2.0) * pz[2] * x3 * x6 - pz[3] * x7 -
                       pz[4] * x10 + (1.0 / 2.0) * pz[5] * x1 * x8 + (1.0 / 2.0) * pz[6] * x6 * x8 - pz[7] * x11;
                J[7] = -pz[0] * x14 - pz[1] * x16 + (1.0 / 2.0) * pz[2] * x15 * x3 + (1.0 / 2.0) * pz[3] * x13 * x3 -
                       pz[4] * x17 - pz[5] * x18 + (1.0 / 2.0) * pz[6] * x15 * x8 + (1.0 / 2.0) * pz[7] * x13 * x8;
                J[8] = -pz[0] * x20 - pz[1] * x21 - pz[2] * x23 - pz[3] * x24 + (1.0 / 2.0) * pz[4] * x1 * x13 +
                       (1.0 / 2.0) * pz[5] * x1 * x15 + (1.0 / 2.0) * pz[6] * x15 * x6 + (1.0 / 2.0) * pz[7] * x13 * x6;
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
                //	- Result: 9*ADD + 9*ASSIGNMENT + 27*MUL
                //	- Subexpressions: 32*ADD + 44*DIV + 132*MUL + 42*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * x3;
                T x5 = x1 * x4;
                T x6 = x0 + 0.5;
                T x7 = x4 * x6;
                T x8 = x2 + 0.5;
                T x9 = (1.0 / 2.0) * x8;
                T x10 = x6 * x9;
                T x11 = x1 * x9;
                T x12 = -py[0] * x5 - py[1] * x7 - py[2] * x10 - py[3] * x11 + (1.0 / 2.0) * py[4] * x1 * x3 +
                        (1.0 / 2.0) * py[5] * x3 * x6 + (1.0 / 2.0) * py[6] * x6 * x8 + (1.0 / 2.0) * py[7] * x1 * x8;
                T x13 = (1.0 / 2.0) * z;
                T x14 = 0.5 - x13;
                T x15 = (1.0 / 2.0) * x14;
                T x16 = x1 * x15;
                T x17 = x15 * x6;
                T x18 = x13 + 0.5;
                T x19 = (1.0 / 2.0) * x18;
                T x20 = x1 * x19;
                T x21 = x19 * x6;
                T x22 = -pz[0] * x16 - pz[1] * x17 + (1.0 / 2.0) * pz[2] * x14 * x6 + (1.0 / 2.0) * pz[3] * x1 * x14 -
                        pz[4] * x20 - pz[5] * x21 + (1.0 / 2.0) * pz[6] * x18 * x6 + (1.0 / 2.0) * pz[7] * x1 * x18;
                T x23 = x12 * x22;
                T x24 = -py[0] * x16 - py[1] * x17 + (1.0 / 2.0) * py[2] * x14 * x6 + (1.0 / 2.0) * py[3] * x1 * x14 -
                        py[4] * x20 - py[5] * x21 + (1.0 / 2.0) * py[6] * x18 * x6 + (1.0 / 2.0) * py[7] * x1 * x18;
                T x25 = -pz[0] * x5 - pz[1] * x7 - pz[2] * x10 - pz[3] * x11 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                        (1.0 / 2.0) * pz[5] * x3 * x6 + (1.0 / 2.0) * pz[6] * x6 * x8 + (1.0 / 2.0) * pz[7] * x1 * x8;
                T x26 = x15 * x3;
                T x27 = x15 * x8;
                T x28 = x19 * x3;
                T x29 = x19 * x8;
                T x30 = -pz[0] * x26 + (1.0 / 2.0) * pz[1] * x14 * x3 + (1.0 / 2.0) * pz[2] * x14 * x8 - pz[3] * x27 -
                        pz[4] * x28 + (1.0 / 2.0) * pz[5] * x18 * x3 + (1.0 / 2.0) * pz[6] * x18 * x8 - pz[7] * x29;
                T x31 = -px[0] * x5 - px[1] * x7 - px[2] * x10 - px[3] * x11 + (1.0 / 2.0) * px[4] * x1 * x3 +
                        (1.0 / 2.0) * px[5] * x3 * x6 + (1.0 / 2.0) * px[6] * x6 * x8 + (1.0 / 2.0) * px[7] * x1 * x8;
                T x32 = x24 * x31;
                T x33 = -py[0] * x26 + (1.0 / 2.0) * py[1] * x14 * x3 + (1.0 / 2.0) * py[2] * x14 * x8 - py[3] * x27 -
                        py[4] * x28 + (1.0 / 2.0) * py[5] * x18 * x3 + (1.0 / 2.0) * py[6] * x18 * x8 - py[7] * x29;
                T x34 = -px[0] * x16 - px[1] * x17 + (1.0 / 2.0) * px[2] * x14 * x6 + (1.0 / 2.0) * px[3] * x1 * x14 -
                        px[4] * x20 - px[5] * x21 + (1.0 / 2.0) * px[6] * x18 * x6 + (1.0 / 2.0) * px[7] * x1 * x18;
                T x35 = x25 * x34;
                T x36 = -px[0] * x26 + (1.0 / 2.0) * px[1] * x14 * x3 + (1.0 / 2.0) * px[2] * x14 * x8 - px[3] * x27 -
                        px[4] * x28 + (1.0 / 2.0) * px[5] * x18 * x3 + (1.0 / 2.0) * px[6] * x18 * x8 - px[7] * x29;
                T x37 = 1.0 / (x12 * x30 * x34 + x22 * x31 * x33 - x23 * x36 + x24 * x25 * x36 - x30 * x32 - x33 * x35);
                J_inv[0] = x37 * (-x23 + x24 * x25);
                J_inv[1] = x37 * (x22 * x31 - x35);
                J_inv[2] = x37 * (x12 * x34 - x32);
                J_inv[3] = x37 * (x12 * x30 - x25 * x33);
                J_inv[4] = x37 * (x25 * x36 - x30 * x31);
                J_inv[5] = x37 * (-x12 * x36 + x31 * x33);
                J_inv[6] = x37 * (x22 * x33 - x24 * x30);
                J_inv[7] = x37 * (-x22 * x36 + x30 * x34);
                J_inv[8] = x37 * (x24 * x36 - x33 * x34);
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
                //	- Result: 3*ADD + 3*ASSIGNMENT + 24*MUL
                //	- Subexpressions: 3*ADD + 3*DIV + 12*MUL + 3*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * z;
                T x5 = 0.5 - x4;
                T x6 = x3 * x5;
                T x7 = x1 * x6;
                T x8 = x0 + 0.5;
                T x9 = x6 * x8;
                T x10 = x2 + 0.5;
                T x11 = x10 * x5;
                T x12 = x11 * x8;
                T x13 = x1 * x11;
                T x14 = x4 + 0.5;
                T x15 = x14 * x3;
                T x16 = x1 * x15;
                T x17 = x15 * x8;
                T x18 = x10 * x14;
                T x19 = x18 * x8;
                T x20 = x1 * x18;
                tx = px[0] * x7 + px[1] * x9 + px[2] * x12 + px[3] * x13 + px[4] * x16 + px[5] * x17 + px[6] * x19 +
                     px[7] * x20;
                ty = py[0] * x7 + py[1] * x9 + py[2] * x12 + py[3] * x13 + py[4] * x16 + py[5] * x17 + py[6] * x19 +
                     py[7] * x20;
                tz = pz[0] * x7 + pz[1] * x9 + pz[2] * x12 + pz[3] * x13 + pz[4] * x16 + pz[5] * x17 + pz[6] * x19 +
                     pz[7] * x20;
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
                {
                    static constexpr int max_iter = 100;
                    GeoT tol = 1e-16;
                    bool converged = false;
                    GeoT residual_norm;
                    for (int i = 0; i < max_iter; ++i) {
                        // FLOATING POINT OPS!
                        //	- Result: 13*ADD + 4*ASSIGNMENT + 27*MUL + 3*POW
                        //	- Subexpressions: 32*ADD + 44*DIV + 171*MUL + 66*SUB
                        T x0 = (1.0 / 2.0) * x;
                        T x1 = 0.5 - x0;
                        T x2 = (1.0 / 2.0) * y;
                        T x3 = 0.5 - x2;
                        T x4 = (1.0 / 2.0) * x3;
                        T x5 = x1 * x4;
                        T x6 = x0 + 0.5;
                        T x7 = x4 * x6;
                        T x8 = x2 + 0.5;
                        T x9 = (1.0 / 2.0) * x8;
                        T x10 = x6 * x9;
                        T x11 = x1 * x9;
                        T x12 = -px[0] * x5 - px[1] * x7 - px[2] * x10 - px[3] * x11 + (1.0 / 2.0) * px[4] * x1 * x3 +
                                (1.0 / 2.0) * px[5] * x3 * x6 + (1.0 / 2.0) * px[6] * x6 * x8 +
                                (1.0 / 2.0) * px[7] * x1 * x8;
                        T x13 = (1.0 / 2.0) * z;
                        T x14 = 0.5 - x13;
                        T x15 = (1.0 / 2.0) * x14;
                        T x16 = x1 * x15;
                        T x17 = x15 * x6;
                        T x18 = x13 + 0.5;
                        T x19 = (1.0 / 2.0) * x18;
                        T x20 = x1 * x19;
                        T x21 = x19 * x6;
                        T x22 = -py[0] * x16 - py[1] * x17 + (1.0 / 2.0) * py[2] * x14 * x6 +
                                (1.0 / 2.0) * py[3] * x1 * x14 - py[4] * x20 - py[5] * x21 +
                                (1.0 / 2.0) * py[6] * x18 * x6 + (1.0 / 2.0) * py[7] * x1 * x18;
                        T x23 = x12 * x22;
                        T x24 = -px[0] * x16 - px[1] * x17 + (1.0 / 2.0) * px[2] * x14 * x6 +
                                (1.0 / 2.0) * px[3] * x1 * x14 - px[4] * x20 - px[5] * x21 +
                                (1.0 / 2.0) * px[6] * x18 * x6 + (1.0 / 2.0) * px[7] * x1 * x18;
                        T x25 = -py[0] * x5 - py[1] * x7 - py[2] * x10 - py[3] * x11 + (1.0 / 2.0) * py[4] * x1 * x3 +
                                (1.0 / 2.0) * py[5] * x3 * x6 + (1.0 / 2.0) * py[6] * x6 * x8 +
                                (1.0 / 2.0) * py[7] * x1 * x8;
                        T x26 = x15 * x3;
                        T x27 = x15 * x8;
                        T x28 = x19 * x3;
                        T x29 = x19 * x8;
                        T x30 = -pz[0] * x26 + (1.0 / 2.0) * pz[1] * x14 * x3 + (1.0 / 2.0) * pz[2] * x14 * x8 -
                                pz[3] * x27 - pz[4] * x28 + (1.0 / 2.0) * pz[5] * x18 * x3 +
                                (1.0 / 2.0) * pz[6] * x18 * x8 - pz[7] * x29;
                        T x31 = -py[0] * x26 + (1.0 / 2.0) * py[1] * x14 * x3 + (1.0 / 2.0) * py[2] * x14 * x8 -
                                py[3] * x27 - py[4] * x28 + (1.0 / 2.0) * py[5] * x18 * x3 +
                                (1.0 / 2.0) * py[6] * x18 * x8 - py[7] * x29;
                        T x32 = -pz[0] * x5 - pz[1] * x7 - pz[2] * x10 - pz[3] * x11 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                                (1.0 / 2.0) * pz[5] * x3 * x6 + (1.0 / 2.0) * pz[6] * x6 * x8 +
                                (1.0 / 2.0) * pz[7] * x1 * x8;
                        T x33 = x24 * x32;
                        T x34 = -px[0] * x26 + (1.0 / 2.0) * px[1] * x14 * x3 + (1.0 / 2.0) * px[2] * x14 * x8 -
                                px[3] * x27 - px[4] * x28 + (1.0 / 2.0) * px[5] * x18 * x3 +
                                (1.0 / 2.0) * px[6] * x18 * x8 - px[7] * x29;
                        T x35 = -pz[0] * x16 - pz[1] * x17 + (1.0 / 2.0) * pz[2] * x14 * x6 +
                                (1.0 / 2.0) * pz[3] * x1 * x14 - pz[4] * x20 - pz[5] * x21 +
                                (1.0 / 2.0) * pz[6] * x18 * x6 + (1.0 / 2.0) * pz[7] * x1 * x18;
                        T x36 = x25 * x35;
                        T x37 = 1.0 / (x12 * x31 * x35 + x22 * x32 * x34 - x23 * x30 + x24 * x25 * x30 - x31 * x33 -
                                       x34 * x36);
                        T x38 = x14 * x3;
                        T x39 = x1 * x38;
                        T x40 = x38 * x6;
                        T x41 = x14 * x8;
                        T x42 = x41 * x6;
                        T x43 = x1 * x41;
                        T x44 = x18 * x3;
                        T x45 = x1 * x44;
                        T x46 = x44 * x6;
                        T x47 = x18 * x8;
                        T x48 = x47 * x6;
                        T x49 = x1 * x47;
                        T x50 = -pz[0] * x39 - pz[1] * x40 - pz[2] * x42 - pz[3] * x43 - pz[4] * x45 - pz[5] * x46 -
                                pz[6] * x48 - pz[7] * x49 + tz;
                        T x51 = x37 * x50;
                        T x52 = -py[0] * x39 - py[1] * x40 - py[2] * x42 - py[3] * x43 - py[4] * x45 - py[5] * x46 -
                                py[6] * x48 - py[7] * x49 + ty;
                        T x53 = x37 * x52;
                        T x54 = -px[0] * x39 - px[1] * x40 - px[2] * x42 - px[3] * x43 - px[4] * x45 - px[5] * x46 -
                                px[6] * x48 - px[7] * x49 + tx;
                        T x55 = x37 * x54;
                        x = tx + x51 * (-x23 + x24 * x25) + x53 * (x12 * x35 - x33) + x55 * (x22 * x32 - x36);
                        y = ty + x51 * (x12 * x31 - x25 * x34) + x53 * (-x12 * x30 + x32 * x34) +
                            x55 * (x25 * x30 - x31 * x32);
                        z = tz + x51 * (x22 * x34 - x24 * x31) + x53 * (x24 * x30 - x34 * x35) +
                            x55 * (-x22 * x30 + x31 * x35);
                        residual_norm = pow(x50, 2) + pow(x52, 2) + pow(x54, 2);

                        // FIXME not SIMD friendly
                        if (residual_norm < tol) {
                            converged = true;
                            break;
                        }
                    }

                    assert(converged);
                }
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
                //	- Result: 24*ADD + 24*ASSIGNMENT + 69*MUL
                //	- Subexpressions: 35*ADD + 47*DIV + 177*MUL + 5*NEG + 53*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * x3;
                T x5 = x1 * x4;
                T x6 = x0 + 0.5;
                T x7 = x4 * x6;
                T x8 = x2 + 0.5;
                T x9 = (1.0 / 2.0) * x8;
                T x10 = x6 * x9;
                T x11 = x1 * x9;
                T x12 = -px[0] * x5 - px[1] * x7 - px[2] * x10 - px[3] * x11 + (1.0 / 2.0) * px[4] * x1 * x3 +
                        (1.0 / 2.0) * px[5] * x3 * x6 + (1.0 / 2.0) * px[6] * x6 * x8 + (1.0 / 2.0) * px[7] * x1 * x8;
                T x13 = (1.0 / 2.0) * z;
                T x14 = 0.5 - x13;
                T x15 = (1.0 / 2.0) * x14;
                T x16 = x1 * x15;
                T x17 = x15 * x6;
                T x18 = x13 + 0.5;
                T x19 = (1.0 / 2.0) * x18;
                T x20 = x1 * x19;
                T x21 = x19 * x6;
                T x22 = -pz[0] * x16 - pz[1] * x17 + (1.0 / 2.0) * pz[2] * x14 * x6 + (1.0 / 2.0) * pz[3] * x1 * x14 -
                        pz[4] * x20 - pz[5] * x21 + (1.0 / 2.0) * pz[6] * x18 * x6 + (1.0 / 2.0) * pz[7] * x1 * x18;
                T x23 = -px[0] * x16 - px[1] * x17 + (1.0 / 2.0) * px[2] * x14 * x6 + (1.0 / 2.0) * px[3] * x1 * x14 -
                        px[4] * x20 - px[5] * x21 + (1.0 / 2.0) * px[6] * x18 * x6 + (1.0 / 2.0) * px[7] * x1 * x18;
                T x24 = -pz[0] * x5 - pz[1] * x7 - pz[2] * x10 - pz[3] * x11 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                        (1.0 / 2.0) * pz[5] * x3 * x6 + (1.0 / 2.0) * pz[6] * x6 * x8 + (1.0 / 2.0) * pz[7] * x1 * x8;
                T x25 = x23 * x24;
                T x26 = x12 * x22 - x25;
                T x27 = (1.0 / 4.0) * x;
                T x28 = x27 - 0.25;
                T x29 = x15 * x3;
                T x30 = x15 * x8;
                T x31 = x19 * x3;
                T x32 = x19 * x8;
                T x33 = -pz[0] * x29 + (1.0 / 2.0) * pz[1] * x14 * x3 + (1.0 / 2.0) * pz[2] * x14 * x8 - pz[3] * x30 -
                        pz[4] * x31 + (1.0 / 2.0) * pz[5] * x18 * x3 + (1.0 / 2.0) * pz[6] * x18 * x8 - pz[7] * x32;
                T x34 = -py[0] * x16 - py[1] * x17 + (1.0 / 2.0) * py[2] * x14 * x6 + (1.0 / 2.0) * py[3] * x1 * x14 -
                        py[4] * x20 - py[5] * x21 + (1.0 / 2.0) * py[6] * x18 * x6 + (1.0 / 2.0) * py[7] * x1 * x18;
                T x35 = x12 * x34;
                T x36 = -py[0] * x29 + (1.0 / 2.0) * py[1] * x14 * x3 + (1.0 / 2.0) * py[2] * x14 * x8 - py[3] * x30 -
                        py[4] * x31 + (1.0 / 2.0) * py[5] * x18 * x3 + (1.0 / 2.0) * py[6] * x18 * x8 - py[7] * x32;
                T x37 = -px[0] * x29 + (1.0 / 2.0) * px[1] * x14 * x3 + (1.0 / 2.0) * px[2] * x14 * x8 - px[3] * x30 -
                        px[4] * x31 + (1.0 / 2.0) * px[5] * x18 * x3 + (1.0 / 2.0) * px[6] * x18 * x8 - px[7] * x32;
                T x38 = -py[0] * x5 - py[1] * x7 - py[2] * x10 - py[3] * x11 + (1.0 / 2.0) * py[4] * x1 * x3 +
                        (1.0 / 2.0) * py[5] * x3 * x6 + (1.0 / 2.0) * py[6] * x6 * x8 + (1.0 / 2.0) * py[7] * x1 * x8;
                T x39 = x22 * x38;
                T x40 = 1.0 / (x12 * x22 * x36 + x23 * x33 * x38 + x24 * x34 * x37 - x25 * x36 - x33 * x35 - x37 * x39);
                T x41 = x14 * x40;
                T x42 = x28 * x41;
                T x43 = x24 * x34 - x39;
                T x44 = (1.0 / 4.0) * y;
                T x45 = x44 - 0.25;
                T x46 = x41 * x45;
                T x47 = x23 * x38 - x35;
                T x48 = x40 * x5;
                T x49 = -x12 * x33 + x24 * x37;
                T x50 = -x24 * x36 + x33 * x38;
                T x51 = x12 * x36 - x37 * x38;
                T x52 = -x22 * x37 + x23 * x33;
                T x53 = x22 * x36 - x33 * x34;
                T x54 = -x23 * x36 + x34 * x37;
                T x55 = -x41 * x45;
                T x56 = x27 + 0.25;
                T x57 = -x41 * x56;
                T x58 = x40 * x7;
                T x59 = x41 * x56;
                T x60 = x44 + 0.25;
                T x61 = x41 * x60;
                T x62 = x10 * x40;
                T x63 = x47 * x62;
                T x64 = x51 * x62;
                T x65 = x54 * x62;
                T x66 = -x28;
                T x67 = x41 * x66;
                T x68 = -x41 * x60;
                T x69 = x11 * x40;
                T x70 = x3 * x40;
                T x71 = x66 * x70;
                T x72 = (1.0 / 4.0) * z + 0.25;
                T x73 = -x70 * x72;
                T x74 = x20 * x40;
                T x75 = x56 * x70;
                T x76 = x70 * x72;
                T x77 = x21 * x40;
                T x78 = x18 * x40;
                T x79 = x56 * x78;
                T x80 = x60 * x78;
                T x81 = x1 * x40;
                T x82 = x60 * x81;
                T x83 = x72 * x81;
                T x84 = x32 * x40;
                gx[0] = x26 * x42 + x43 * x46 - x47 * x48;
                gy[0] = x42 * x49 + x46 * x50 - x48 * x51;
                gz[0] = x42 * x52 + x46 * x53 - x48 * x54;
                gx[1] = x26 * x57 + x43 * x55 - x47 * x58;
                gy[1] = x49 * x57 + x50 * x55 - x51 * x58;
                gz[1] = x52 * x57 + x53 * x55 - x54 * x58;
                gx[2] = x26 * x59 + x43 * x61 - x63;
                gy[2] = x49 * x59 + x50 * x61 - x64;
                gz[2] = x52 * x59 + x53 * x61 - x65;
                gx[3] = x26 * x67 + x43 * x68 - x47 * x69;
                gy[3] = x49 * x67 + x50 * x68 - x51 * x69;
                gz[3] = x52 * x67 + x53 * x68 - x54 * x69;
                gx[4] = -x26 * x74 + x43 * x73 + x47 * x71;
                gy[4] = -x49 * x74 + x50 * x73 + x51 * x71;
                gz[4] = -x52 * x74 + x53 * x73 + x54 * x71;
                gx[5] = -x26 * x77 + x43 * x76 + x47 * x75;
                gy[5] = -x49 * x77 + x50 * x76 + x51 * x75;
                gz[5] = -x52 * x77 + x53 * x76 + x54 * x75;
                gx[6] = x26 * x79 + x43 * x80 + x63;
                gy[6] = x49 * x79 + x50 * x80 + x64;
                gz[6] = x52 * x79 + x53 * x80 + x65;
                gx[7] = x26 * x83 - x43 * x84 + x47 * x82;
                gy[7] = x49 * x83 - x50 * x84 + x51 * x82;
                gz[7] = x52 * x83 - x53 * x84 + x54 * x82;
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
                //	- Result: 24*ADD + 33*ASSIGNMENT + 77*MUL
                //	- Subexpressions: 35*ADD + 51*DIV + 177*MUL + 5*NEG + 53*SUB
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
                T x13 = (1.0 / 2.0) * x6;
                T x14 = (1.0 / 2.0) * x9;
                T x15 = (1.0 / 2.0) * x11;
                T x16 = (1.0 / 2.0) * x12;
                T x17 = -pz[0] * x13 + (1.0 / 2.0) * pz[1] * x3 * x5 + (1.0 / 2.0) * pz[2] * x5 * x8 - pz[3] * x14 -
                        pz[4] * x15 + (1.0 / 2.0) * pz[5] * x10 * x3 + (1.0 / 2.0) * pz[6] * x10 * x8 - pz[7] * x16;
                T x18 = (1.0 / 2.0) * x3;
                T x19 = x1 * x18;
                T x20 = x18 * x7;
                T x21 = (1.0 / 2.0) * x8;
                T x22 = x21 * x7;
                T x23 = x1 * x21;
                T x24 = -px[0] * x19 - px[1] * x20 - px[2] * x22 - px[3] * x23 + (1.0 / 2.0) * px[4] * x1 * x3 +
                        (1.0 / 2.0) * px[5] * x3 * x7 + (1.0 / 2.0) * px[6] * x7 * x8 + (1.0 / 2.0) * px[7] * x1 * x8;
                T x25 = (1.0 / 2.0) * x5;
                T x26 = x1 * x25;
                T x27 = x25 * x7;
                T x28 = (1.0 / 2.0) * x10;
                T x29 = x1 * x28;
                T x30 = x28 * x7;
                T x31 = -py[0] * x26 - py[1] * x27 + (1.0 / 2.0) * py[2] * x5 * x7 + (1.0 / 2.0) * py[3] * x1 * x5 -
                        py[4] * x29 - py[5] * x30 + (1.0 / 2.0) * py[6] * x10 * x7 + (1.0 / 2.0) * py[7] * x1 * x10;
                T x32 = x24 * x31;
                T x33 = -py[0] * x13 + (1.0 / 2.0) * py[1] * x3 * x5 + (1.0 / 2.0) * py[2] * x5 * x8 - py[3] * x14 -
                        py[4] * x15 + (1.0 / 2.0) * py[5] * x10 * x3 + (1.0 / 2.0) * py[6] * x10 * x8 - py[7] * x16;
                T x34 = -px[0] * x26 - px[1] * x27 + (1.0 / 2.0) * px[2] * x5 * x7 + (1.0 / 2.0) * px[3] * x1 * x5 -
                        px[4] * x29 - px[5] * x30 + (1.0 / 2.0) * px[6] * x10 * x7 + (1.0 / 2.0) * px[7] * x1 * x10;
                T x35 = -pz[0] * x19 - pz[1] * x20 - pz[2] * x22 - pz[3] * x23 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                        (1.0 / 2.0) * pz[5] * x3 * x7 + (1.0 / 2.0) * pz[6] * x7 * x8 + (1.0 / 2.0) * pz[7] * x1 * x8;
                T x36 = x34 * x35;
                T x37 = -px[0] * x13 + (1.0 / 2.0) * px[1] * x3 * x5 + (1.0 / 2.0) * px[2] * x5 * x8 - px[3] * x14 -
                        px[4] * x15 + (1.0 / 2.0) * px[5] * x10 * x3 + (1.0 / 2.0) * px[6] * x10 * x8 - px[7] * x16;
                T x38 = -py[0] * x19 - py[1] * x20 - py[2] * x22 - py[3] * x23 + (1.0 / 2.0) * py[4] * x1 * x3 +
                        (1.0 / 2.0) * py[5] * x3 * x7 + (1.0 / 2.0) * py[6] * x7 * x8 + (1.0 / 2.0) * py[7] * x1 * x8;
                T x39 = -pz[0] * x26 - pz[1] * x27 + (1.0 / 2.0) * pz[2] * x5 * x7 + (1.0 / 2.0) * pz[3] * x1 * x5 -
                        pz[4] * x29 - pz[5] * x30 + (1.0 / 2.0) * pz[6] * x10 * x7 + (1.0 / 2.0) * pz[7] * x1 * x10;
                T x40 = x38 * x39;
                T x41 = -x17 * x32 + x17 * x34 * x38 + x24 * x33 * x39 + x31 * x35 * x37 - x33 * x36 - x37 * x40;
                T x42 = x24 * x39 - x36;
                T x43 = (1.0 / 4.0) * x;
                T x44 = x43 - 0.25;
                T x45 = 1.0 / x41;
                T x46 = x45 * x5;
                T x47 = x44 * x46;
                T x48 = x31 * x35 - x40;
                T x49 = (1.0 / 4.0) * y;
                T x50 = x49 - 0.25;
                T x51 = x46 * x50;
                T x52 = -x32 + x34 * x38;
                T x53 = x19 * x45;
                T x54 = -x17 * x24 + x35 * x37;
                T x55 = x17 * x38 - x33 * x35;
                T x56 = x24 * x33 - x37 * x38;
                T x57 = x17 * x34 - x37 * x39;
                T x58 = -x17 * x31 + x33 * x39;
                T x59 = x31 * x37 - x33 * x34;
                T x60 = -x46 * x50;
                T x61 = x43 + 0.25;
                T x62 = -x46 * x61;
                T x63 = x20 * x45;
                T x64 = x46 * x61;
                T x65 = x49 + 0.25;
                T x66 = x46 * x65;
                T x67 = x22 * x45;
                T x68 = x52 * x67;
                T x69 = x56 * x67;
                T x70 = x59 * x67;
                T x71 = -x44;
                T x72 = x46 * x71;
                T x73 = -x46 * x65;
                T x74 = x23 * x45;
                T x75 = x3 * x45;
                T x76 = x71 * x75;
                T x77 = (1.0 / 4.0) * z + 0.25;
                T x78 = -x75 * x77;
                T x79 = x29 * x45;
                T x80 = x61 * x75;
                T x81 = x75 * x77;
                T x82 = x30 * x45;
                T x83 = x10 * x45;
                T x84 = x61 * x83;
                T x85 = x65 * x83;
                T x86 = x1 * x45;
                T x87 = x65 * x86;
                T x88 = x77 * x86;
                T x89 = x16 * x45;
                f[0] = x1 * x6;
                f[1] = x6 * x7;
                f[2] = x7 * x9;
                f[3] = x1 * x9;
                f[4] = x1 * x11;
                f[5] = x11 * x7;
                f[6] = x12 * x7;
                f[7] = x1 * x12;
                measure_value = x41;
                gx[0] = x42 * x47 + x48 * x51 - x52 * x53;
                gy[0] = x47 * x54 + x51 * x55 - x53 * x56;
                gz[0] = x47 * x57 + x51 * x58 - x53 * x59;
                gx[1] = x42 * x62 + x48 * x60 - x52 * x63;
                gy[1] = x54 * x62 + x55 * x60 - x56 * x63;
                gz[1] = x57 * x62 + x58 * x60 - x59 * x63;
                gx[2] = x42 * x64 + x48 * x66 - x68;
                gy[2] = x54 * x64 + x55 * x66 - x69;
                gz[2] = x57 * x64 + x58 * x66 - x70;
                gx[3] = x42 * x72 + x48 * x73 - x52 * x74;
                gy[3] = x54 * x72 + x55 * x73 - x56 * x74;
                gz[3] = x57 * x72 + x58 * x73 - x59 * x74;
                gx[4] = -x42 * x79 + x48 * x78 + x52 * x76;
                gy[4] = -x54 * x79 + x55 * x78 + x56 * x76;
                gz[4] = -x57 * x79 + x58 * x78 + x59 * x76;
                gx[5] = -x42 * x82 + x48 * x81 + x52 * x80;
                gy[5] = -x54 * x82 + x55 * x81 + x56 * x80;
                gz[5] = -x57 * x82 + x58 * x81 + x59 * x80;
                gx[6] = x42 * x84 + x48 * x85 + x68;
                gy[6] = x54 * x84 + x55 * x85 + x69;
                gz[6] = x57 * x84 + x58 * x85 + x70;
                gx[7] = x42 * x88 - x48 * x89 + x52 * x87;
                gy[7] = x54 * x88 - x55 * x89 + x56 * x87;
                gz[7] = x57 * x88 - x58 * x89 + x59 * x87;
            }
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_FE_Hex8_3_IMPL_hpp
