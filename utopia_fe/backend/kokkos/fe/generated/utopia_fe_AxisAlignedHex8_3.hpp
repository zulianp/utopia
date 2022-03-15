#ifndef UTOPIA_TPL_FE_AxisAlignedHex8_3_IMPL_hpp
#define UTOPIA_TPL_FE_AxisAlignedHex8_3_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"

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
		template<typename T, typename GeoT = T>
		class AxisAlignedHex8 {
		public:
			static constexpr int Dim = 3;
			static constexpr int NNodes = 8;
			static constexpr int Order = 1;

			using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

			UTOPIA_FUNCTION static constexpr const char* class_name() { return "AxisAlignedHex8"; }

			UTOPIA_INLINE_FUNCTION static constexpr int dim() 
			{
				return Dim;
			}

			UTOPIA_INLINE_FUNCTION static constexpr int n_nodes() 
			{
				return NNodes;
			}

			UTOPIA_INLINE_FUNCTION static constexpr int order() 
			{
				return Order;
			}

			UTOPIA_FUNCTION static constexpr Result measure(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				// Input quadrature point
				const T x,
				const T y,
				const T z)
			{
				T measure_value;
				//FLOATING POINT OPS!
//	- Result: 3*ADD + ASSIGNMENT + 4*MUL
//	- Subexpressions: 0
measure_value = (-px[0] + px[6])*(-py[0] + py[6])*(-pz[0] + pz[6]);
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
				GeoT *UTOPIA_RESTRICT J)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
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
				GeoT *UTOPIA_RESTRICT J_inv)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 3*ADD + 9*ASSIGNMENT + 6*MUL + 3*POW
//	- Subexpressions: 0
J_inv[0] = -1/(px[0] - px[6]);
J_inv[1] = 0;
J_inv[2] = 0;
J_inv[3] = 0;
J_inv[4] = -1/(py[0] - py[6]);
J_inv[5] = 0;
J_inv[6] = 0;
J_inv[7] = 0;
J_inv[8] = -1/(pz[0] - pz[6]);
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
				GeoT &tz)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 6*ADD + 3*ASSIGNMENT + 6*MUL
//	- Subexpressions: 0
tx = px[0] - x*(px[0] - px[6]);
ty = py[0] - y*(py[0] - py[6]);
tz = pz[0] - z*(pz[0] - pz[6]);
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
				GeoT &z)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 6*ADD + 3*ASSIGNMENT + 9*MUL + 3*POW
//	- Subexpressions: 0
x = (px[0] - tx)/(px[0] - px[6]);
y = (py[0] - ty)/(py[0] - py[6]);
z = (pz[0] - tz)/(pz[0] - pz[6]);
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
				Result *UTOPIA_RESTRICT gz)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 24*ASSIGNMENT + 12*MUL
//	- Subexpressions: 3*DIV + 18*MUL + 6*SUB
T x0 = y - 1.0;
T x1 = 1.0/(px[0] - px[6]);
T x2 = z - 1.0;
T x3 = x1*x2;
T x4 = x0*x3;
T x5 = x - 1.0;
T x6 = 1.0/(py[0] - py[6]);
T x7 = x2*x6;
T x8 = x5*x7;
T x9 = 1.0/(pz[0] - pz[6]);
T x10 = x0*x9;
T x11 = x10*x5;
T x12 = x*x7;
T x13 = x*x10;
T x14 = x3*y;
T x15 = x9*y;
T x16 = x*x15;
T x17 = x15*x5;
T x18 = x1*z;
T x19 = x0*x18;
T x20 = x6*z;
T x21 = x20*x5;
T x22 = x*x20;
T x23 = x18*y;
gx[0] = x4;
gy[0] = x8;
gz[0] = x11;
gx[1] = -x4;
gy[1] = -x12;
gz[1] = -x13;
gx[2] = x14;
gy[2] = x12;
gz[2] = x16;
gx[3] = -x14;
gy[3] = -x8;
gz[3] = -x17;
gx[4] = -x19;
gy[4] = -x21;
gz[4] = -x11;
gx[5] = x19;
gy[5] = x22;
gz[5] = x13;
gx[6] = -x23;
gy[6] = -x22;
gz[6] = -x16;
gx[7] = x23;
gy[7] = x21;
gz[7] = x17;
			}

			UTOPIA_FUNCTION static void value(
				const T x,
				const T y,
				const T z,
				Result *UTOPIA_RESTRICT f
				)
			{
				using namespace utopia::device;
			    // Automatically generated
				//FLOATING POINT OPS!
//	- Result: 8*ASSIGNMENT + 8*MUL
//	- Subexpressions: 4*MUL + 3*SUB
T x0 = 1.0 - x;
T x1 = 1.0 - y;
T x2 = 1.0 - z;
T x3 = x1*x2;
T x4 = x2*y;
T x5 = x1*z;
T x6 = y*z;
f[0] = x0*x3;
f[1] = x*x3;
f[2] = x*x4;
f[3] = x0*x4;
f[4] = x0*x5;
f[5] = x*x5;
f[6] = x*x6;
f[7] = x0*x6;
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
				T &measure_value)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 33*ASSIGNMENT + 21*MUL
//	- Subexpressions: 3*DIV + 22*MUL + 3*NEG + 6*SUB
T x0 = x - 1.0;
T x1 = -x0;
T x2 = y - 1.0;
T x3 = -x2;
T x4 = z - 1.0;
T x5 = -x4;
T x6 = x3*x5;
T x7 = x5*y;
T x8 = x3*z;
T x9 = y*z;
T x10 = px[0] - px[6];
T x11 = py[0] - py[6];
T x12 = pz[0] - pz[6];
T x13 = 1.0/x10;
T x14 = x13*x4;
T x15 = x14*x2;
T x16 = 1.0/x11;
T x17 = x16*x4;
T x18 = x0*x17;
T x19 = 1.0/x12;
T x20 = x19*x2;
T x21 = x0*x20;
T x22 = x*x17;
T x23 = x*x20;
T x24 = x14*y;
T x25 = x19*y;
T x26 = x*x25;
T x27 = x0*x25;
T x28 = x13*x2*z;
T x29 = x16*z;
T x30 = x0*x29;
T x31 = x*x29;
T x32 = x13*x9;
f[0] = x1*x6;
f[1] = x*x6;
f[2] = x*x7;
f[3] = x1*x7;
f[4] = x1*x8;
f[5] = x*x8;
f[6] = x*x9;
f[7] = x1*x9;
measure_value = -x10*x11*x12;
gx[0] = x15;
gy[0] = x18;
gz[0] = x21;
gx[1] = -x15;
gy[1] = -x22;
gz[1] = -x23;
gx[2] = x24;
gy[2] = x22;
gz[2] = x26;
gx[3] = -x24;
gy[3] = -x18;
gz[3] = -x27;
gx[4] = -x28;
gy[4] = -x30;
gz[4] = -x21;
gx[5] = x28;
gy[5] = x31;
gz[5] = x23;
gx[6] = -x32;
gy[6] = -x31;
gz[6] = -x26;
gx[7] = x32;
gy[7] = x30;
gz[7] = x27;
			}


		};
	}
}

#endif // UTOPIA_TPL_FE_AxisAlignedHex8_3_IMPL_hpp
