#ifndef UTOPIA_TPL_FE_Tri3_2_IMPL_hpp
#define UTOPIA_TPL_FE_Tri3_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

// #include "utopia_fe_Tri3.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Tri3 for dimension 2
         */
        template <typename T, typename GeoT>
        class Tri3 {
        public:
            static constexpr int Dim = 2;
            static constexpr int NNodes = 3;
            static constexpr int Order = 1;

            using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Tri3"; }

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

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                // FLOATING POINT OPS!
                //	- Result: 5*ADD + ASSIGNMENT + 6*MUL
                //	- Subexpressions: 0
                measure_value = (-px[0] + px[1]) * (-py[0] + py[2]) - (-px[0] + px[2]) * (-py[0] + py[1]);
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

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                // FLOATING POINT OPS!
                //	- Result: 4*ADD + 4*ASSIGNMENT + 4*MUL
                //	- Subexpressions: 0
                J[0] = -px[0] + px[1];
                J[1] = -py[0] + py[1];
                J[2] = -px[0] + px[2];
                J[3] = -py[0] + py[2];
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

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                // FLOATING POINT OPS!
                //	- Result: 4*ASSIGNMENT + 4*MUL
                //	- Subexpressions: DIV + 2*MUL + 5*SUB
                T x0 = py[0] - py[2];
                T x1 = px[0] - px[1];
                T x2 = px[0] - px[2];
                T x3 = py[0] - py[1];
                T x4 = 1.0 / (x0 * x1 - x2 * x3);
                J_inv[0] = -x0 * x4;
                J_inv[1] = x3 * x4;
                J_inv[2] = x2 * x4;
                J_inv[3] = -x1 * x4;
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
                //	- Result: 2*ADD + 2*ASSIGNMENT + 6*MUL
                //	- Subexpressions: 2*SUB
                T x0 = -x - y + 1;
                tx = px[0] * x0 + px[1] * x + px[2] * y;
                ty = py[0] * x0 + py[1] * x + py[2] * y;
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
                //	- Result: 2*ADD + 2*ASSIGNMENT + 9*MUL
                //	- Subexpressions: 2*ADD + DIV + 8*MUL + 4*SUB
                T x0 = px[0] * py[2];
                T x1 = -px[0] * ty + py[0] * tx;
                T x2 = px[0] * py[1] - px[1] * py[0];
                T x3 = 1.0 / (px[1] * py[2] + px[2] * py[0] - px[2] * py[1] - x0 + x2);
                x = x3 * (px[2] * py[0] - px[2] * ty + py[2] * tx - x0 - x1);
                y = x3 * (px[1] * ty - py[1] * tx + x1 + x2);
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

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                // FLOATING POINT OPS!
                //	- Result: 2*ADD + 6*ASSIGNMENT + 8*MUL
                //	- Subexpressions: 2*ADD + DIV + 2*MUL + 2*NEG + 3*SUB
                T x0 = -py[2];
                T x1 = px[0] - px[1];
                T x2 = py[0] + x0;
                T x3 = -px[2];
                T x4 = px[0] + x3;
                T x5 = py[0] - py[1];
                T x6 = 1.0 / (x1 * x2 - x4 * x5);
                gx[0] = x6 * (py[1] + x0);
                gy[0] = x6 * (-px[1] - x3);
                gx[1] = -x2 * x6;
                gy[1] = x4 * x6;
                gx[2] = x5 * x6;
                gy[2] = -x1 * x6;
            }

            UTOPIA_FUNCTION static void value(const T x, const T y, Result *UTOPIA_RESTRICT f) {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: ADD + 3*ASSIGNMENT + 2*MUL
                //	- Subexpressions: 0
                f[0] = -x - y + 1;
                f[1] = x;
                f[2] = y;
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
                //	- Result: 4*ADD + 10*ASSIGNMENT + 12*MUL
                //	- Subexpressions: 2*ADD + DIV + 2*MUL + 4*NEG + 3*SUB
                T x0 = px[0] - px[1];
                T x1 = -x0;
                T x2 = -py[2];
                T x3 = py[0] + x2;
                T x4 = -x3;
                T x5 = -px[2];
                T x6 = px[0] + x5;
                T x7 = py[0] - py[1];
                T x8 = 1.0 / (x0 * x3 - x6 * x7);
                f[0] = -x - y + 1;
                f[1] = x;
                f[2] = y;
                measure_value = x1 * x4 - x6 * x7;
                gx[0] = x8 * (py[1] + x2);
                gy[0] = x8 * (-px[1] - x5);
                gx[1] = x4 * x8;
                gy[1] = x6 * x8;
                gx[2] = x7 * x8;
                gy[2] = x1 * x8;
            }
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_FE_Tri3_2_IMPL_hpp
