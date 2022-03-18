#ifndef UTOPIA_TPL_MATERIAL_LaplaceOperator_Quad4_2_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_LaplaceOperator_Quad4_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Quad4_2.hpp"
#include "utopia_material_LaplaceOperator.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of LaplaceOperator for symmetric element pair trial=test=Quad4
         */
        template <typename T, typename GeoT>
        class LaplaceOperator<Quad4<T, GeoT>> {
        public:
            using ElemT = Quad4<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "LaplaceOperator<Quad4>"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    // TODO
                }

                // TODO
            };

            LaplaceOperator(const Params &params = Params()) {
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
                //	- Result: 4*ADD + 16*ADDAUGMENTEDASSIGNMENT + 12*MUL + 8*POW
                //	- Subexpressions: 20*ADD + 16*DIV + 53*MUL + 9*NEG + POW + 24*SUB
                T x0 = x - 1.0;
                T x1 = px[0] * x0;
                T x2 = x + 1.0;
                T x3 = px[2] * x2;
                T x4 = px[1] * x2;
                T x5 = px[3] * x0;
                T x6 = x1 + x3 - x4 - x5;
                T x7 = x0 * x6;
                T x8 = y - 1.0;
                T x9 = py[0] * x0;
                T x10 = py[2] * x2;
                T x11 = py[1] * x2;
                T x12 = py[3] * x0;
                T x13 = x10 - x11 - x12 + x9;
                T x14 = -x13 * x8;
                T x15 = x14 + x7;
                T x16 = py[0] * x8;
                T x17 = y + 1.0;
                T x18 = py[2] * x17;
                T x19 = py[1] * x8;
                T x20 = py[3] * x17;
                T x21 = x16 + x18 - x19 - x20;
                T x22 = px[0] * x8;
                T x23 = px[2] * x17;
                T x24 = px[1] * x8;
                T x25 = px[3] * x17;
                T x26 = x22 + x23 - x24 - x25;
                T x27 = pow(-x13 * x26 + x21 * x6, -2);
                T x28 = x0 * x26;
                T x29 = -x21 * x8;
                T x30 = -x28 - x29;
                T x31 = weight * (-((1.0 / 4.0) * x1 + (1.0 / 4.0) * x3 - 1.0 / 4.0 * x4 - 1.0 / 4.0 * x5) *
                                      ((1.0 / 4.0) * x16 + (1.0 / 4.0) * x18 - 1.0 / 4.0 * x19 - 1.0 / 4.0 * x20) +
                                  ((1.0 / 4.0) * x10 - 1.0 / 4.0 * x11 - 1.0 / 4.0 * x12 + (1.0 / 4.0) * x9) *
                                      ((1.0 / 4.0) * x22 + (1.0 / 4.0) * x23 - 1.0 / 4.0 * x24 - 1.0 / 4.0 * x25));
                T x32 = x2 * x6;
                T x33 = -x14 - x32;
                T x34 = x15 * x27;
                T x35 = x2 * x26;
                T x36 = x29 + x35;
                T x37 = x27 * x30;
                T x38 = x31 * (x33 * x34 + x36 * x37);
                T x39 = -x13 * x17;
                T x40 = x32 + x39;
                T x41 = -x17 * x21;
                T x42 = -x35 - x41;
                T x43 = x31 * (x34 * x40 + x37 * x42);
                T x44 = -x39 - x7;
                T x45 = x28 + x41;
                T x46 = x31 * (x34 * x44 + x37 * x45);
                T x47 = x27 * x33;
                T x48 = x27 * x36;
                T x49 = x31 * (x40 * x47 + x42 * x48);
                T x50 = x31 * (x44 * x47 + x45 * x48);
                T x51 = x31 * (x27 * x40 * x44 + x27 * x42 * x45);
                H[0] += x31 * (pow(x15, 2) * x27 + x27 * pow(x30, 2));
                H[1] += x38;
                H[2] += x43;
                H[3] += x46;
                H[4] += x38;
                H[5] += x31 * (x27 * pow(x33, 2) + x27 * pow(x36, 2));
                H[6] += x49;
                H[7] += x50;
                H[8] += x43;
                H[9] += x49;
                H[10] += x31 * (x27 * pow(x40, 2) + x27 * pow(x42, 2));
                H[11] += x51;
                H[12] += x46;
                H[13] += x50;
                H[14] += x51;
                H[15] += x31 * (x27 * pow(x44, 2) + x27 * pow(x45, 2));
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
                //	- Result: 4*ADD + 4*ADDAUGMENTEDASSIGNMENT + 12*MUL
                //	- Subexpressions: 20*ADD + 17*DIV + 45*MUL + 8*NEG + 24*SUB
                T x0 = x - 1.0;
                T x1 = px[0] * x0;
                T x2 = x + 1.0;
                T x3 = px[2] * x2;
                T x4 = px[1] * x2;
                T x5 = px[3] * x0;
                T x6 = x1 + x3 - x4 - x5;
                T x7 = y - 1.0;
                T x8 = py[0] * x7;
                T x9 = y + 1.0;
                T x10 = py[2] * x9;
                T x11 = py[1] * x7;
                T x12 = py[3] * x9;
                T x13 = x10 - x11 - x12 + x8;
                T x14 = px[0] * x7;
                T x15 = px[2] * x9;
                T x16 = px[1] * x7;
                T x17 = px[3] * x9;
                T x18 = x14 + x15 - x16 - x17;
                T x19 = py[0] * x0;
                T x20 = py[2] * x2;
                T x21 = py[1] * x2;
                T x22 = py[3] * x0;
                T x23 = x19 + x20 - x21 - x22;
                T x24 = 1.0 / (x13 * x6 - x18 * x23);
                T x25 = x0 * x6;
                T x26 = -x23 * x7;
                T x27 = x24 * (x25 + x26);
                T x28 = x2 * x6;
                T x29 = -x26 - x28;
                T x30 = u[1] * x24;
                T x31 = -x23 * x9;
                T x32 = x28 + x31;
                T x33 = u[2] * x24;
                T x34 = -x25 - x31;
                T x35 = u[3] * x24;
                T x36 = u[0] * x27 + x29 * x30 + x32 * x33 + x34 * x35;
                T x37 = x0 * x18;
                T x38 = -x13 * x7;
                T x39 = x24 * (-x37 - x38);
                T x40 = x18 * x2;
                T x41 = x38 + x40;
                T x42 = -x13 * x9;
                T x43 = -x40 - x42;
                T x44 = x37 + x42;
                T x45 = u[0] * x39 + x30 * x41 + x33 * x43 + x35 * x44;
                T x46 = 4 * weight *
                        (-((1.0 / 4.0) * x1 + (1.0 / 4.0) * x3 - 1.0 / 4.0 * x4 - 1.0 / 4.0 * x5) *
                             ((1.0 / 4.0) * x10 - 1.0 / 4.0 * x11 - 1.0 / 4.0 * x12 + (1.0 / 4.0) * x8) +
                         ((1.0 / 4.0) * x14 + (1.0 / 4.0) * x15 - 1.0 / 4.0 * x16 - 1.0 / 4.0 * x17) *
                             ((1.0 / 4.0) * x19 + (1.0 / 4.0) * x20 - 1.0 / 4.0 * x21 - 1.0 / 4.0 * x22));
                T x47 = x24 * x36;
                T x48 = x24 * x45;
                Hx[0] += x46 * (x27 * x36 + x39 * x45);
                Hx[1] += x46 * (x29 * x47 + x41 * x48);
                Hx[2] += x46 * (x32 * x47 + x43 * x48);
                Hx[3] += x46 * (x34 * x47 + x44 * x48);
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
                //	- Result: 4*ADD + 4*ADDAUGMENTEDASSIGNMENT + 12*MUL
                //	- Subexpressions: 20*ADD + 17*DIV + 45*MUL + 8*NEG + 24*SUB
                T x0 = x - 1.0;
                T x1 = px[0] * x0;
                T x2 = x + 1.0;
                T x3 = px[2] * x2;
                T x4 = px[1] * x2;
                T x5 = px[3] * x0;
                T x6 = x1 + x3 - x4 - x5;
                T x7 = y - 1.0;
                T x8 = py[0] * x7;
                T x9 = y + 1.0;
                T x10 = py[2] * x9;
                T x11 = py[1] * x7;
                T x12 = py[3] * x9;
                T x13 = x10 - x11 - x12 + x8;
                T x14 = px[0] * x7;
                T x15 = px[2] * x9;
                T x16 = px[1] * x7;
                T x17 = px[3] * x9;
                T x18 = x14 + x15 - x16 - x17;
                T x19 = py[0] * x0;
                T x20 = py[2] * x2;
                T x21 = py[1] * x2;
                T x22 = py[3] * x0;
                T x23 = x19 + x20 - x21 - x22;
                T x24 = 1.0 / (x13 * x6 - x18 * x23);
                T x25 = x0 * x6;
                T x26 = -x23 * x7;
                T x27 = x24 * (x25 + x26);
                T x28 = x2 * x6;
                T x29 = -x26 - x28;
                T x30 = u[1] * x24;
                T x31 = -x23 * x9;
                T x32 = x28 + x31;
                T x33 = u[2] * x24;
                T x34 = -x25 - x31;
                T x35 = u[3] * x24;
                T x36 = u[0] * x27 + x29 * x30 + x32 * x33 + x34 * x35;
                T x37 = x0 * x18;
                T x38 = -x13 * x7;
                T x39 = x24 * (-x37 - x38);
                T x40 = x18 * x2;
                T x41 = x38 + x40;
                T x42 = -x13 * x9;
                T x43 = -x40 - x42;
                T x44 = x37 + x42;
                T x45 = u[0] * x39 + x30 * x41 + x33 * x43 + x35 * x44;
                T x46 = 4 * weight *
                        (-((1.0 / 4.0) * x1 + (1.0 / 4.0) * x3 - 1.0 / 4.0 * x4 - 1.0 / 4.0 * x5) *
                             ((1.0 / 4.0) * x10 - 1.0 / 4.0 * x11 - 1.0 / 4.0 * x12 + (1.0 / 4.0) * x8) +
                         ((1.0 / 4.0) * x14 + (1.0 / 4.0) * x15 - 1.0 / 4.0 * x16 - 1.0 / 4.0 * x17) *
                             ((1.0 / 4.0) * x19 + (1.0 / 4.0) * x20 - 1.0 / 4.0 * x21 - 1.0 / 4.0 * x22));
                T x47 = x24 * x36;
                T x48 = x24 * x45;
                g[0] += x46 * (x27 * x36 + x39 * x45);
                g[1] += x46 * (x29 * x47 + x41 * x48);
                g[2] += x46 * (x32 * x47 + x43 * x48);
                g[3] += x46 * (x34 * x47 + x44 * x48);
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
                //	- Result: 16*ADD + ADDAUGMENTEDASSIGNMENT + 35*MUL + 2*POW
                //	- Subexpressions: 6*ADD + DIV + 30*MUL + 4*NEG + 11*SUB
                T x0 = x - 1.0;
                T x1 = px[0] * x0;
                T x2 = x + 1.0;
                T x3 = px[1] * x2;
                T x4 = px[2] * x2;
                T x5 = px[3] * x0;
                T x6 = y - 1.0;
                T x7 = py[0] * x6;
                T x8 = py[1] * x6;
                T x9 = y + 1.0;
                T x10 = py[2] * x9;
                T x11 = py[3] * x9;
                T x12 = px[0] * x6;
                T x13 = px[1] * x6;
                T x14 = px[2] * x9;
                T x15 = px[3] * x9;
                T x16 = py[0] * x0;
                T x17 = py[1] * x2;
                T x18 = py[2] * x2;
                T x19 = py[3] * x0;
                T x20 = x1 - x3 + x4 - x5;
                T x21 = x0 * x20;
                T x22 = x16 - x17 + x18 - x19;
                T x23 = -x22 * x6;
                T x24 = x10 - x11 + x7 - x8;
                T x25 = x12 - x13 + x14 - x15;
                T x26 = 1.0 / (x20 * x24 - x22 * x25);
                T x27 = u[0] * x26;
                T x28 = x2 * x20;
                T x29 = u[1] * x26;
                T x30 = -x22 * x9;
                T x31 = u[2] * x26;
                T x32 = u[3] * x26;
                T x33 = x0 * x25;
                T x34 = -x24 * x6;
                T x35 = x2 * x25;
                T x36 = -x24 * x9;
                e += weight *
                     (-((1.0 / 4.0) * x1 - 1.0 / 4.0 * x3 + (1.0 / 4.0) * x4 - 1.0 / 4.0 * x5) *
                          ((1.0 / 4.0) * x10 - 1.0 / 4.0 * x11 + (1.0 / 4.0) * x7 - 1.0 / 4.0 * x8) +
                      ((1.0 / 4.0) * x12 - 1.0 / 4.0 * x13 + (1.0 / 4.0) * x14 - 1.0 / 4.0 * x15) *
                          ((1.0 / 4.0) * x16 - 1.0 / 4.0 * x17 + (1.0 / 4.0) * x18 - 1.0 / 4.0 * x19)) *
                     (pow(x27 * (x21 + x23) + x29 * (-x23 - x28) + x31 * (x28 + x30) + x32 * (-x21 - x30), 2) +
                      pow(x27 * (-x33 - x34) + x29 * (x34 + x35) + x31 * (-x35 - x36) + x32 * (x33 + x36), 2));
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
                //	- Result: 9*ADD + 21*ADDAUGMENTEDASSIGNMENT + 25*MUL + 10*POW
                //	- Subexpressions: 26*ADD + 17*DIV + 69*MUL + 9*NEG + POW + 24*SUB
                T x0 = x - 1.0;
                T x1 = px[0] * x0;
                T x2 = x + 1.0;
                T x3 = px[2] * x2;
                T x4 = px[1] * x2;
                T x5 = px[3] * x0;
                T x6 = x1 + x3 - x4 - x5;
                T x7 = x0 * x6;
                T x8 = y - 1.0;
                T x9 = py[0] * x0;
                T x10 = py[2] * x2;
                T x11 = py[1] * x2;
                T x12 = py[3] * x0;
                T x13 = x10 - x11 - x12 + x9;
                T x14 = -x13 * x8;
                T x15 = x14 + x7;
                T x16 = py[0] * x8;
                T x17 = y + 1.0;
                T x18 = py[2] * x17;
                T x19 = py[1] * x8;
                T x20 = py[3] * x17;
                T x21 = x16 + x18 - x19 - x20;
                T x22 = px[0] * x8;
                T x23 = px[2] * x17;
                T x24 = px[1] * x8;
                T x25 = px[3] * x17;
                T x26 = x22 + x23 - x24 - x25;
                T x27 = -x13 * x26 + x21 * x6;
                T x28 = pow(x27, -2);
                T x29 = x0 * x26;
                T x30 = -x21 * x8;
                T x31 = -x29 - x30;
                T x32 = weight * (-((1.0 / 4.0) * x1 + (1.0 / 4.0) * x3 - 1.0 / 4.0 * x4 - 1.0 / 4.0 * x5) *
                                      ((1.0 / 4.0) * x16 + (1.0 / 4.0) * x18 - 1.0 / 4.0 * x19 - 1.0 / 4.0 * x20) +
                                  ((1.0 / 4.0) * x10 - 1.0 / 4.0 * x11 - 1.0 / 4.0 * x12 + (1.0 / 4.0) * x9) *
                                      ((1.0 / 4.0) * x22 + (1.0 / 4.0) * x23 - 1.0 / 4.0 * x24 - 1.0 / 4.0 * x25));
                T x33 = x2 * x6;
                T x34 = -x14 - x33;
                T x35 = x15 * x28;
                T x36 = x2 * x26;
                T x37 = x30 + x36;
                T x38 = x28 * x31;
                T x39 = x32 * (x34 * x35 + x37 * x38);
                T x40 = -x13 * x17;
                T x41 = x33 + x40;
                T x42 = -x17 * x21;
                T x43 = -x36 - x42;
                T x44 = x32 * (x35 * x41 + x38 * x43);
                T x45 = -x40 - x7;
                T x46 = x29 + x42;
                T x47 = x32 * (x35 * x45 + x38 * x46);
                T x48 = x28 * x34;
                T x49 = x28 * x37;
                T x50 = x32 * (x41 * x48 + x43 * x49);
                T x51 = x32 * (x45 * x48 + x46 * x49);
                T x52 = x32 * (x28 * x41 * x45 + x28 * x43 * x46);
                T x53 = 1.0 / x27;
                T x54 = x15 * x53;
                T x55 = u[1] * x53;
                T x56 = u[2] * x53;
                T x57 = u[3] * x53;
                T x58 = u[0] * x54 + x34 * x55 + x41 * x56 + x45 * x57;
                T x59 = x31 * x53;
                T x60 = u[0] * x59 + x37 * x55 + x43 * x56 + x46 * x57;
                T x61 = 4 * x32;
                T x62 = x53 * x58;
                T x63 = x53 * x60;
                H[0] += x32 * (pow(x15, 2) * x28 + x28 * pow(x31, 2));
                H[1] += x39;
                H[2] += x44;
                H[3] += x47;
                H[4] += x39;
                H[5] += x32 * (x28 * pow(x34, 2) + x28 * pow(x37, 2));
                H[6] += x50;
                H[7] += x51;
                H[8] += x44;
                H[9] += x50;
                H[10] += x32 * (x28 * pow(x41, 2) + x28 * pow(x43, 2));
                H[11] += x52;
                H[12] += x47;
                H[13] += x51;
                H[14] += x52;
                H[15] += x32 * (x28 * pow(x45, 2) + x28 * pow(x46, 2));
                g[0] += x61 * (x54 * x58 + x59 * x60);
                g[1] += x61 * (x34 * x62 + x37 * x63);
                g[2] += x61 * (x41 * x62 + x43 * x63);
                g[3] += x61 * (x45 * x62 + x46 * x63);
                e += x32 * (pow(x58, 2) + pow(x60, 2));
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using LaplaceOperatorQuad4 = utopia::kokkos::AutoKernel<
            FE,
            utopia::kernels::LaplaceOperator<utopia::kernels::Quad4<typename FE::Scalar, typename FE::Scalar>>,
            2>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_LaplaceOperator_Quad4_2_IMPL_hpp
