#ifndef UTOPIA_TPL_FE_AxisAlignedQuad4_2_IMPL_hpp
#define UTOPIA_TPL_FE_AxisAlignedQuad4_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

// #include "utopia_fe_AxisAlignedQuad4.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of AxisAlignedQuad4 for dimension 2
         */
        template <typename T, typename GeoT>
        class AxisAlignedQuad4 {
        public:
            static constexpr int Dim = 2;
            static constexpr int NNodes = 4;
            static constexpr int Order = 1;

            using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "AxisAlignedQuad4"; }

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
                //	- Result: 2*ADD + ASSIGNMENT + 3*MUL
                //	- Subexpressions: 0
                measure_value = (-px[0] + px[2]) * (-py[0] + py[2]);
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
                //	- Result: 2*ADD + 4*ASSIGNMENT + 2*MUL
                //	- Subexpressions: 0
                J[0] = -px[0] + px[2];
                J[1] = 0;
                J[2] = 0;
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
                //	- Result: 2*ADD + 4*ASSIGNMENT + 4*MUL + 2*POW
                //	- Subexpressions: 0
                J_inv[0] = -1 / (px[0] - px[2]);
                J_inv[1] = 0;
                J_inv[2] = 0;
                J_inv[3] = -1 / (py[0] - py[2]);
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
                //	- Result: 4*ADD + 2*ASSIGNMENT + 4*MUL
                //	- Subexpressions: 0
                tx = px[0] - x * (px[0] - px[2]);
                ty = py[0] - y * (py[0] - py[2]);
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
                //	- Result: 4*ADD + 2*ASSIGNMENT + 6*MUL + 2*POW
                //	- Subexpressions: 0
                x = (px[0] - tx) / (px[0] - px[2]);
                y = (py[0] - ty) / (py[0] - py[2]);
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
                //	- Result: 8*ASSIGNMENT + 8*MUL
                //	- Subexpressions: 2*ADD + 2*DIV + 4*SUB
                T x0 = y - 1;
                T x1 = 0.25 / (px[0] - px[2]);
                T x2 = x - 1;
                T x3 = 0.25 / (py[0] - py[2]);
                T x4 = x + 1;
                T x5 = y + 1;
                gx[0] = -x0 * x1;
                gy[0] = -x2 * x3;
                gx[1] = x0 * x1;
                gy[1] = x3 * x4;
                gx[2] = -x1 * x5;
                gy[2] = -x3 * x4;
                gx[3] = x1 * x5;
                gy[3] = x2 * x3;
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
                //	- Result: 13*ASSIGNMENT + 13*MUL
                //	- Subexpressions: 4*ADD + 4*DIV + 6*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x0 + 0.5;
                T x5 = x2 + 0.5;
                T x6 = px[0] - px[2];
                T x7 = py[0] - py[2];
                T x8 = y - 1;
                T x9 = 0.25 / x6;
                T x10 = x - 1;
                T x11 = 0.25 / x7;
                T x12 = x + 1;
                T x13 = y + 1;
                f[0] = x1 * x3;
                f[1] = x3 * x4;
                f[2] = x4 * x5;
                f[3] = x1 * x5;
                measure_value = x6 * x7;
                gx[0] = -x8 * x9;
                gy[0] = -x10 * x11;
                gx[1] = x8 * x9;
                gy[1] = x11 * x12;
                gx[2] = -x13 * x9;
                gy[2] = -x11 * x12;
                gx[3] = x13 * x9;
                gy[3] = x10 * x11;
            }
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_FE_AxisAlignedQuad4_2_IMPL_hpp
