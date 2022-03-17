#ifndef UTOPIA_TPL_MATERIAL_Mass_Pentatope5_4_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_Pentatope5_4_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Pentatope5_4.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Mass for symmetric element pair trial=test=Pentatope5
         */
        template <typename T, typename GeoT>
        class Mass<Pentatope5<T, GeoT>> {
        public:
            using ElemT = Pentatope5<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Mass<Pentatope5>"; }

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
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 25*ADDAUGMENTEDASSIGNMENT + 5*MUL + 5*POW
                //	- Subexpressions: 13*ADD + 56*MUL + NEG + 30*SUB
                T x0 = -t - x - y - z + 1;
                T x1 = px[0] - px[2];
                T x2 = -x1;
                T x3 = -py[0] + py[4];
                T x4 = -pz[0] + pz[3];
                T x5 = x3 * x4;
                T x6 = -px[0] + px[3];
                T x7 = -pz[0] + pz[4];
                T x8 = -py[0] + py[2];
                T x9 = x7 * x8;
                T x10 = -px[0] + px[4];
                T x11 = -py[0] + py[3];
                T x12 = -pz[0] + pz[2];
                T x13 = x11 * x12;
                T x14 = x11 * x7;
                T x15 = x12 * x3;
                T x16 = x4 * x8;
                T x17 = -pt[0] + pt[2];
                T x18 = x10 * x17;
                T x19 = -pt[0] + pt[3];
                T x20 = x19 * x2;
                T x21 = -pt[0] + pt[4];
                T x22 = x21 * x6;
                T x23 = x17 * x6;
                T x24 = x10 * x19;
                T x25 = weight *
                        ((-pt[0] + pt[1]) * (x1 * x14 + x10 * x13 - x10 * x16 - x15 * x6 + x2 * x5 + x6 * x9) +
                         (-px[0] + px[1]) * (-x13 * x21 + x14 * x17 + x15 * x19 + x16 * x21 - x17 * x5 - x19 * x9) +
                         (-py[0] + py[1]) * (x1 * x21 * x4 + x12 * x22 - x12 * x24 + x18 * x4 + x20 * x7 - x23 * x7) +
                         (-pz[0] + pz[1]) * (-x11 * x18 + x11 * x2 * x21 - x20 * x3 - x22 * x8 + x23 * x3 + x24 * x8));
                T x26 = x0 * x25;
                T x27 = x * x26;
                T x28 = x26 * y;
                T x29 = x26 * z;
                T x30 = t * x26;
                T x31 = x * x25;
                T x32 = x31 * y;
                T x33 = x31 * z;
                T x34 = t * x31;
                T x35 = x25 * y;
                T x36 = x35 * z;
                T x37 = t * x35;
                T x38 = t * x25 * z;
                H[0] += pow(x0, 2) * x25;
                H[1] += x27;
                H[2] += x28;
                H[3] += x29;
                H[4] += x30;
                H[5] += x27;
                H[6] += pow(x, 2) * x25;
                H[7] += x32;
                H[8] += x33;
                H[9] += x34;
                H[10] += x28;
                H[11] += x32;
                H[12] += x25 * pow(y, 2);
                H[13] += x36;
                H[14] += x37;
                H[15] += x29;
                H[16] += x33;
                H[17] += x36;
                H[18] += x25 * pow(z, 2);
                H[19] += x38;
                H[20] += x30;
                H[21] += x34;
                H[22] += x37;
                H[23] += x38;
                H[24] += pow(t, 2) * x25;
            }

            UTOPIA_FUNCTION void apply(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T *UTOPIA_RESTRICT Hx) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 5*ADDAUGMENTEDASSIGNMENT + 5*MUL
                //	- Subexpressions: 17*ADD + 48*MUL + NEG + 30*SUB
                T x0 = -t - x - y - z + 1;
                T x1 = px[0] - px[2];
                T x2 = -x1;
                T x3 = -py[0] + py[4];
                T x4 = -pz[0] + pz[3];
                T x5 = x3 * x4;
                T x6 = -px[0] + px[3];
                T x7 = -pz[0] + pz[4];
                T x8 = -py[0] + py[2];
                T x9 = x7 * x8;
                T x10 = -px[0] + px[4];
                T x11 = -py[0] + py[3];
                T x12 = -pz[0] + pz[2];
                T x13 = x11 * x12;
                T x14 = x11 * x7;
                T x15 = x12 * x3;
                T x16 = x4 * x8;
                T x17 = -pt[0] + pt[2];
                T x18 = x10 * x17;
                T x19 = -pt[0] + pt[3];
                T x20 = x19 * x2;
                T x21 = -pt[0] + pt[4];
                T x22 = x21 * x6;
                T x23 = x17 * x6;
                T x24 = x10 * x19;
                T x25 = weight *
                        ((-pt[0] + pt[1]) * (x1 * x14 + x10 * x13 - x10 * x16 - x15 * x6 + x2 * x5 + x6 * x9) +
                         (-px[0] + px[1]) * (-x13 * x21 + x14 * x17 + x15 * x19 + x16 * x21 - x17 * x5 - x19 * x9) +
                         (-py[0] + py[1]) * (x1 * x21 * x4 + x12 * x22 - x12 * x24 + x18 * x4 + x20 * x7 - x23 * x7) +
                         (-pz[0] + pz[1]) * (-x11 * x18 + x11 * x2 * x21 - x20 * x3 - x22 * x8 + x23 * x3 + x24 * x8)) *
                        (t * u[4] + u[0] * x0 + u[1] * x + u[2] * y + u[3] * z);
                Hx[0] += x0 * x25;
                Hx[1] += x * x25;
                Hx[2] += x25 * y;
                Hx[3] += x25 * z;
                Hx[4] += t * x25;
            }

            UTOPIA_FUNCTION void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T *UTOPIA_RESTRICT g) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 5*ADDAUGMENTEDASSIGNMENT + 5*MUL
                //	- Subexpressions: 17*ADD + 48*MUL + NEG + 30*SUB
                T x0 = -t - x - y - z + 1;
                T x1 = px[0] - px[2];
                T x2 = -x1;
                T x3 = -py[0] + py[4];
                T x4 = -pz[0] + pz[3];
                T x5 = x3 * x4;
                T x6 = -px[0] + px[3];
                T x7 = -pz[0] + pz[4];
                T x8 = -py[0] + py[2];
                T x9 = x7 * x8;
                T x10 = -px[0] + px[4];
                T x11 = -py[0] + py[3];
                T x12 = -pz[0] + pz[2];
                T x13 = x11 * x12;
                T x14 = x11 * x7;
                T x15 = x12 * x3;
                T x16 = x4 * x8;
                T x17 = -pt[0] + pt[2];
                T x18 = x10 * x17;
                T x19 = -pt[0] + pt[3];
                T x20 = x19 * x2;
                T x21 = -pt[0] + pt[4];
                T x22 = x21 * x6;
                T x23 = x17 * x6;
                T x24 = x10 * x19;
                T x25 = weight *
                        ((-pt[0] + pt[1]) * (x1 * x14 + x10 * x13 - x10 * x16 - x15 * x6 + x2 * x5 + x6 * x9) +
                         (-px[0] + px[1]) * (-x13 * x21 + x14 * x17 + x15 * x19 + x16 * x21 - x17 * x5 - x19 * x9) +
                         (-py[0] + py[1]) * (x1 * x21 * x4 + x12 * x22 - x12 * x24 + x18 * x4 + x20 * x7 - x23 * x7) +
                         (-pz[0] + pz[1]) * (-x11 * x18 + x11 * x2 * x21 - x20 * x3 - x22 * x8 + x23 * x3 + x24 * x8)) *
                        (t * u[4] + u[0] * x0 + u[1] * x + u[2] * y + u[3] * z);
                g[0] += x0 * x25;
                g[1] += x * x25;
                g[2] += x25 * y;
                g[3] += x25 * z;
                g[4] += t * x25;
            }

            UTOPIA_FUNCTION void value(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T &e) const {
                using namespace utopia::device;
                // Automatically generated
                // TODO
            }

            UTOPIA_FUNCTION void eval(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T &e,
                T *UTOPIA_RESTRICT H,
                T *UTOPIA_RESTRICT g) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 30*ADDAUGMENTEDASSIGNMENT + 10*MUL + 5*POW
                //	- Subexpressions: 17*ADD + 61*MUL + NEG + 30*SUB
                T x0 = -t - x - y - z + 1;
                T x1 = px[0] - px[2];
                T x2 = -x1;
                T x3 = -py[0] + py[4];
                T x4 = -pz[0] + pz[3];
                T x5 = x3 * x4;
                T x6 = -px[0] + px[3];
                T x7 = -pz[0] + pz[4];
                T x8 = -py[0] + py[2];
                T x9 = x7 * x8;
                T x10 = -px[0] + px[4];
                T x11 = -py[0] + py[3];
                T x12 = -pz[0] + pz[2];
                T x13 = x11 * x12;
                T x14 = x11 * x7;
                T x15 = x12 * x3;
                T x16 = x4 * x8;
                T x17 = -pt[0] + pt[2];
                T x18 = x10 * x17;
                T x19 = -pt[0] + pt[3];
                T x20 = x19 * x2;
                T x21 = -pt[0] + pt[4];
                T x22 = x21 * x6;
                T x23 = x17 * x6;
                T x24 = x10 * x19;
                T x25 = weight *
                        ((-pt[0] + pt[1]) * (x1 * x14 + x10 * x13 - x10 * x16 - x15 * x6 + x2 * x5 + x6 * x9) +
                         (-px[0] + px[1]) * (-x13 * x21 + x14 * x17 + x15 * x19 + x16 * x21 - x17 * x5 - x19 * x9) +
                         (-py[0] + py[1]) * (x1 * x21 * x4 + x12 * x22 - x12 * x24 + x18 * x4 + x20 * x7 - x23 * x7) +
                         (-pz[0] + pz[1]) * (-x11 * x18 + x11 * x2 * x21 - x20 * x3 - x22 * x8 + x23 * x3 + x24 * x8));
                T x26 = x0 * x25;
                T x27 = x * x26;
                T x28 = x26 * y;
                T x29 = x26 * z;
                T x30 = t * x26;
                T x31 = x * x25;
                T x32 = x31 * y;
                T x33 = x31 * z;
                T x34 = t * x31;
                T x35 = x25 * y;
                T x36 = x35 * z;
                T x37 = t * x35;
                T x38 = x25 * z;
                T x39 = t * x38;
                T x40 = t * u[4] + u[0] * x0 + u[1] * x + u[2] * y + u[3] * z;
                H[0] += pow(x0, 2) * x25;
                H[1] += x27;
                H[2] += x28;
                H[3] += x29;
                H[4] += x30;
                H[5] += x27;
                H[6] += pow(x, 2) * x25;
                H[7] += x32;
                H[8] += x33;
                H[9] += x34;
                H[10] += x28;
                H[11] += x32;
                H[12] += x25 * pow(y, 2);
                H[13] += x36;
                H[14] += x37;
                H[15] += x29;
                H[16] += x33;
                H[17] += x36;
                H[18] += x25 * pow(z, 2);
                H[19] += x39;
                H[20] += x30;
                H[21] += x34;
                H[22] += x37;
                H[23] += x39;
                H[24] += pow(t, 2) * x25;
                g[0] += x26 * x40;
                g[1] += x31 * x40;
                g[2] += x35 * x40;
                g[3] += x38 * x40;
                g[4] += t * x25 * x40;
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using MassPentatope5 = utopia::kokkos::AutoKernel<
            FE,
            utopia::kernels::Mass<utopia::kernels::Pentatope5<typename FE::Scalar, typename FE::Scalar>>,
            4>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_Mass_Pentatope5_4_IMPL_hpp
