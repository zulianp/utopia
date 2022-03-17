#ifndef UTOPIA_TPL_MATERIAL_Mass_Tet4_3_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_Tet4_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Tet4_3.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Mass for symmetric element pair trial=test=Tet4
         */
        template <typename T, typename GeoT>
        class Mass<Tet4<T, GeoT>> {
        public:
            using ElemT = Tet4<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Mass<Tet4>"; }

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
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 16*ADDAUGMENTEDASSIGNMENT + 4*MUL + 4*POW
                //	- Subexpressions: 2*ADD + 22*MUL + 15*SUB
                T x0 = -x - y - z + 1;
                T x1 = -px[0] + px[1];
                T x2 = -py[0] + py[2];
                T x3 = -pz[0] + pz[3];
                T x4 = -px[0] + px[2];
                T x5 = -py[0] + py[3];
                T x6 = -pz[0] + pz[1];
                T x7 = -px[0] + px[3];
                T x8 = -py[0] + py[1];
                T x9 = -pz[0] + pz[2];
                T x10 =
                    weight * (x1 * x2 * x3 - x1 * x5 * x9 - x2 * x6 * x7 - x3 * x4 * x8 + x4 * x5 * x6 + x7 * x8 * x9);
                T x11 = x0 * x10;
                T x12 = x * x11;
                T x13 = x11 * y;
                T x14 = x11 * z;
                T x15 = x * x10;
                T x16 = x15 * y;
                T x17 = x15 * z;
                T x18 = x10 * y * z;
                H[0] += pow(x0, 2) * x10;
                H[1] += x12;
                H[2] += x13;
                H[3] += x14;
                H[4] += x12;
                H[5] += pow(x, 2) * x10;
                H[6] += x16;
                H[7] += x17;
                H[8] += x13;
                H[9] += x16;
                H[10] += x10 * pow(y, 2);
                H[11] += x18;
                H[12] += x14;
                H[13] += x17;
                H[14] += x18;
                H[15] += x10 * pow(z, 2);
            }

            UTOPIA_FUNCTION void apply(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T *UTOPIA_RESTRICT Hx) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 4*ADDAUGMENTEDASSIGNMENT + 4*MUL
                //	- Subexpressions: 5*ADD + 18*MUL + 15*SUB
                T x0 = -x - y - z + 1;
                T x1 = -px[0] + px[1];
                T x2 = -py[0] + py[2];
                T x3 = -pz[0] + pz[3];
                T x4 = -px[0] + px[2];
                T x5 = -py[0] + py[3];
                T x6 = -pz[0] + pz[1];
                T x7 = -px[0] + px[3];
                T x8 = -py[0] + py[1];
                T x9 = -pz[0] + pz[2];
                T x10 = weight * (u[0] * x0 + u[1] * x + u[2] * y + u[3] * z) *
                        (x1 * x2 * x3 - x1 * x5 * x9 - x2 * x6 * x7 - x3 * x4 * x8 + x4 * x5 * x6 + x7 * x8 * x9);
                Hx[0] += x0 * x10;
                Hx[1] += x * x10;
                Hx[2] += x10 * y;
                Hx[3] += x10 * z;
            }

            UTOPIA_FUNCTION void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T *UTOPIA_RESTRICT g) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 4*ADDAUGMENTEDASSIGNMENT + 4*MUL
                //	- Subexpressions: 5*ADD + 18*MUL + 15*SUB
                T x0 = -x - y - z + 1;
                T x1 = -px[0] + px[1];
                T x2 = -py[0] + py[2];
                T x3 = -pz[0] + pz[3];
                T x4 = -px[0] + px[2];
                T x5 = -py[0] + py[3];
                T x6 = -pz[0] + pz[1];
                T x7 = -px[0] + px[3];
                T x8 = -py[0] + py[1];
                T x9 = -pz[0] + pz[2];
                T x10 = weight * (u[0] * x0 + u[1] * x + u[2] * y + u[3] * z) *
                        (x1 * x2 * x3 - x1 * x5 * x9 - x2 * x6 * x7 - x3 * x4 * x8 + x4 * x5 * x6 + x7 * x8 * x9);
                g[0] += x0 * x10;
                g[1] += x * x10;
                g[2] += x10 * y;
                g[3] += x10 * z;
            }

            UTOPIA_FUNCTION void value(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
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
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T &e,
                T *UTOPIA_RESTRICT g,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 20*ADDAUGMENTEDASSIGNMENT + 8*MUL + 4*POW
                //	- Subexpressions: 5*ADD + 26*MUL + 15*SUB
                T x0 = -x - y - z + 1;
                T x1 = -px[0] + px[1];
                T x2 = -py[0] + py[2];
                T x3 = -pz[0] + pz[3];
                T x4 = -px[0] + px[2];
                T x5 = -py[0] + py[3];
                T x6 = -pz[0] + pz[1];
                T x7 = -px[0] + px[3];
                T x8 = -py[0] + py[1];
                T x9 = -pz[0] + pz[2];
                T x10 =
                    weight * (x1 * x2 * x3 - x1 * x5 * x9 - x2 * x6 * x7 - x3 * x4 * x8 + x4 * x5 * x6 + x7 * x8 * x9);
                T x11 = x0 * x10;
                T x12 = x * x11;
                T x13 = x11 * y;
                T x14 = x11 * z;
                T x15 = x * x10;
                T x16 = x15 * y;
                T x17 = x15 * z;
                T x18 = x10 * y;
                T x19 = x18 * z;
                T x20 = u[0] * x0 + u[1] * x + u[2] * y + u[3] * z;
                H[0] += pow(x0, 2) * x10;
                H[1] += x12;
                H[2] += x13;
                H[3] += x14;
                H[4] += x12;
                H[5] += pow(x, 2) * x10;
                H[6] += x16;
                H[7] += x17;
                H[8] += x13;
                H[9] += x16;
                H[10] += x10 * pow(y, 2);
                H[11] += x19;
                H[12] += x14;
                H[13] += x17;
                H[14] += x19;
                H[15] += x10 * pow(z, 2);
                g[0] += x11 * x20;
                g[1] += x15 * x20;
                g[2] += x18 * x20;
                g[3] += x10 * x20 * z;
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using MassTet4 = utopia::kokkos::
            AutoKernel<FE, utopia::kernels::Mass<utopia::kernels::Tet4<typename FE::Scalar, typename FE::Scalar>>, 3>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_Mass_Tet4_3_IMPL_hpp
