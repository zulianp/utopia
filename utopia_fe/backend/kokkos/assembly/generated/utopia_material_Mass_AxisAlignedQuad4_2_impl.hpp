#ifndef UTOPIA_TPL_MATERIAL_Mass_AxisAlignedQuad4_2_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_AxisAlignedQuad4_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_AxisAlignedQuad4_2.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Mass for symmetric element pair trial=test=AxisAlignedQuad4
         */
        template <typename T, typename GeoT>
        class Mass<AxisAlignedQuad4<T, GeoT>> {
        public:
            using ElemT = AxisAlignedQuad4<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Mass<AxisAlignedQuad4>"; }

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
                //	- Subexpressions: 4*ADD + 4*DIV + 14*MUL + 4*POW + 6*SUB
                T x0 = pow(1 - x, 2);
                T x1 = pow(1 - y, 2);
                T x2 = weight * (-px[0] + px[2]) * (-py[0] + py[2]);
                T x3 = (1.0 / 16.0) * x2;
                T x4 = x1 * x3;
                T x5 = (1.0 / 4.0) * x2;
                T x6 = (1.0 / 2.0) * x;
                T x7 = (0.5 - x6) * (x6 + 0.5);
                T x8 = x5 * x7;
                T x9 = x1 * x8;
                T x10 = (1.0 / 2.0) * y;
                T x11 = (0.5 - x10) * (x10 + 0.5);
                T x12 = x11 * x2 * x7;
                T x13 = x11 * x5;
                T x14 = x0 * x13;
                T x15 = pow(x + 1, 2);
                T x16 = x13 * x15;
                T x17 = pow(y + 1, 2);
                T x18 = x17 * x3;
                T x19 = x17 * x8;
                H[0] += x0 * x4;
                H[1] += x9;
                H[2] += x12;
                H[3] += x14;
                H[4] += x9;
                H[5] += x15 * x4;
                H[6] += x16;
                H[7] += x12;
                H[8] += x12;
                H[9] += x16;
                H[10] += x15 * x18;
                H[11] += x19;
                H[12] += x14;
                H[13] += x12;
                H[14] += x19;
                H[15] += x0 * x18;
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
                //	- Subexpressions: 5*ADD + 2*DIV + 11*MUL + 4*SUB
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
                T x10 = weight * (-px[0] + px[2]) * (-py[0] + py[2]) * (u[0] * x4 + u[1] * x6 + u[2] * x8 + u[3] * x9);
                Hx[0] += x10 * x4;
                Hx[1] += x10 * x6;
                Hx[2] += x10 * x8;
                Hx[3] += x10 * x9;
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
                //	- Subexpressions: 5*ADD + 2*DIV + 11*MUL + 4*SUB
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
                T x10 = weight * (-px[0] + px[2]) * (-py[0] + py[2]) * (u[0] * x4 + u[1] * x6 + u[2] * x8 + u[3] * x9);
                g[0] += x10 * x4;
                g[1] += x10 * x6;
                g[2] += x10 * x8;
                g[3] += x10 * x9;
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
                //	- Result: 3*ADD + ADDAUGMENTEDASSIGNMENT + 7*MUL + POW
                //	- Subexpressions: 2*ADD + 2*DIV + 2*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = x0 + 0.5;
                T x5 = x2 + 0.5;
                e += weight * (-px[0] + px[2]) * (-py[0] + py[2]) *
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
                //	- Subexpressions: 7*ADD + 4*DIV + 23*MUL + 4*POW + 6*SUB
                T x0 = pow(1 - x, 2);
                T x1 = pow(1 - y, 2);
                T x2 = weight * (-px[0] + px[2]) * (-py[0] + py[2]);
                T x3 = (1.0 / 16.0) * x2;
                T x4 = x1 * x3;
                T x5 = (1.0 / 2.0) * x;
                T x6 = 0.5 - x5;
                T x7 = x5 + 0.5;
                T x8 = (1.0 / 4.0) * x2;
                T x9 = x6 * x7 * x8;
                T x10 = x1 * x9;
                T x11 = (1.0 / 2.0) * y;
                T x12 = 0.5 - x11;
                T x13 = x12 * x6;
                T x14 = x11 + 0.5;
                T x15 = x14 * x7;
                T x16 = x13 * x15 * x2;
                T x17 = x12 * x14 * x8;
                T x18 = x0 * x17;
                T x19 = pow(x + 1, 2);
                T x20 = x17 * x19;
                T x21 = pow(y + 1, 2);
                T x22 = x21 * x3;
                T x23 = x21 * x9;
                T x24 = x12 * x7;
                T x25 = x14 * x6;
                T x26 = u[0] * x13 + u[1] * x24 + u[2] * x15 + u[3] * x25;
                T x27 = x2 * x26;
                H[0] += x0 * x4;
                H[1] += x10;
                H[2] += x16;
                H[3] += x18;
                H[4] += x10;
                H[5] += x19 * x4;
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
                g[0] += x13 * x27;
                g[1] += x24 * x27;
                g[2] += x15 * x27;
                g[3] += x25 * x27;
                e += x2 * pow(x26, 2);
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FunctionSpace, class FE>
        using MassAxisAlignedQuad4 = utopia::kokkos::AutoKernel<FunctionSpace,
            FE,
            utopia::kernels::Mass<utopia::kernels::AxisAlignedQuad4<typename FE::Scalar, typename FE::Scalar>>,
            2>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_Mass_AxisAlignedQuad4_2_IMPL_hpp
