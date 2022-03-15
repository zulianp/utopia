#ifndef UTOPIA_TPL_FE_Hex8_3_IMPL_hpp
#define UTOPIA_TPL_FE_Hex8_3_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"

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
		template<typename T, typename GeoT = T>
		class Hex8 {
		public:
			static constexpr int Dim = 3;
			static constexpr int NNodes = 8;
			static constexpr int Order = 1;

			using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

			UTOPIA_FUNCTION static constexpr const char* class_name() { return "Hex8"; }

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
//	- Result: ADD + ASSIGNMENT + 6*MUL
//	- Subexpressions: 27*ADD + 120*MUL + 39*SUB
T x0 = x*y;
T x1 = 1.0 - y;
T x2 = x*x1;
T x3 = 1.0 - x;
T x4 = x3*y;
T x5 = x1*x3;
T x6 = -px[0]*x5 - px[1]*x2 - px[2]*x0 - px[3]*x4 + px[4]*x1*x3 + px[5]*x*x1 + px[6]*x*y + px[7]*x3*y;
T x7 = x*z;
T x8 = 1.0 - z;
T x9 = x*x8;
T x10 = x3*z;
T x11 = x3*x8;
T x12 = -py[0]*x11 - py[1]*x9 + py[2]*x*x8 + py[3]*x3*x8 - py[4]*x10 - py[5]*x7 + py[6]*x*z + py[7]*x3*z;
T x13 = y*z;
T x14 = x8*y;
T x15 = x1*z;
T x16 = x1*x8;
T x17 = -pz[0]*x16 + pz[1]*x1*x8 + pz[2]*x8*y - pz[3]*x14 - pz[4]*x15 + pz[5]*x1*z + pz[6]*y*z - pz[7]*x13;
T x18 = -px[0]*x11 - px[1]*x9 + px[2]*x*x8 + px[3]*x3*x8 - px[4]*x10 - px[5]*x7 + px[6]*x*z + px[7]*x3*z;
T x19 = -py[0]*x16 + py[1]*x1*x8 + py[2]*x8*y - py[3]*x14 - py[4]*x15 + py[5]*x1*z + py[6]*y*z - py[7]*x13;
T x20 = -pz[0]*x5 - pz[1]*x2 - pz[2]*x0 - pz[3]*x4 + pz[4]*x1*x3 + pz[5]*x*x1 + pz[6]*x*y + pz[7]*x3*y;
T x21 = -px[0]*x16 + px[1]*x1*x8 + px[2]*x8*y - px[3]*x14 - px[4]*x15 + px[5]*x1*z + px[6]*y*z - px[7]*x13;
T x22 = -py[0]*x5 - py[1]*x2 - py[2]*x0 - py[3]*x4 + py[4]*x1*x3 + py[5]*x*x1 + py[6]*x*y + py[7]*x3*y;
T x23 = -pz[0]*x11 - pz[1]*x9 + pz[2]*x*x8 + pz[3]*x3*x8 - pz[4]*x10 - pz[5]*x7 + pz[6]*x*z + pz[7]*x3*z;
measure_value = -x12*x17*x6 + x12*x20*x21 + x17*x18*x22 - x18*x19*x20 + x19*x23*x6 - x21*x22*x23;
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
//	- Result: 9*ADD + 9*ASSIGNMENT + 72*MUL
//	- Subexpressions: 12*MUL + 3*SUB
T x0 = y*z;
T x1 = 1.0 - z;
T x2 = x1*y;
T x3 = 1.0 - y;
T x4 = x3*z;
T x5 = x1*x3;
T x6 = x*z;
T x7 = x*x1;
T x8 = 1.0 - x;
T x9 = x8*z;
T x10 = x1*x8;
T x11 = x*y;
T x12 = x*x3;
T x13 = x8*y;
T x14 = x3*x8;
J[0] = -px[0]*x5 + px[1]*x1*x3 + px[2]*x1*y - px[3]*x2 - px[4]*x4 + px[5]*x3*z + px[6]*y*z - px[7]*x0;
J[1] = -px[0]*x10 - px[1]*x7 + px[2]*x*x1 + px[3]*x1*x8 - px[4]*x9 - px[5]*x6 + px[6]*x*z + px[7]*x8*z;
J[2] = -px[0]*x14 - px[1]*x12 - px[2]*x11 - px[3]*x13 + px[4]*x3*x8 + px[5]*x*x3 + px[6]*x*y + px[7]*x8*y;
J[3] = -py[0]*x5 + py[1]*x1*x3 + py[2]*x1*y - py[3]*x2 - py[4]*x4 + py[5]*x3*z + py[6]*y*z - py[7]*x0;
J[4] = -py[0]*x10 - py[1]*x7 + py[2]*x*x1 + py[3]*x1*x8 - py[4]*x9 - py[5]*x6 + py[6]*x*z + py[7]*x8*z;
J[5] = -py[0]*x14 - py[1]*x12 - py[2]*x11 - py[3]*x13 + py[4]*x3*x8 + py[5]*x*x3 + py[6]*x*y + py[7]*x8*y;
J[6] = -pz[0]*x5 + pz[1]*x1*x3 + pz[2]*x1*y - pz[3]*x2 - pz[4]*x4 + pz[5]*x3*z + pz[6]*y*z - pz[7]*x0;
J[7] = -pz[0]*x10 - pz[1]*x7 + pz[2]*x*x1 + pz[3]*x1*x8 - pz[4]*x9 - pz[5]*x6 + pz[6]*x*z + pz[7]*x8*z;
J[8] = -pz[0]*x14 - pz[1]*x12 - pz[2]*x11 - pz[3]*x13 + pz[4]*x3*x8 + pz[5]*x*x3 + pz[6]*x*y + pz[7]*x8*y;
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
//	- Result: 9*ADD + 9*ASSIGNMENT + 27*MUL
//	- Subexpressions: 29*ADD + DIV + 132*MUL + 42*SUB
T x0 = x*y;
T x1 = 1.0 - y;
T x2 = x*x1;
T x3 = 1.0 - x;
T x4 = x3*y;
T x5 = x1*x3;
T x6 = -py[0]*x5 - py[1]*x2 - py[2]*x0 - py[3]*x4 + py[4]*x1*x3 + py[5]*x*x1 + py[6]*x*y + py[7]*x3*y;
T x7 = x*z;
T x8 = 1.0 - z;
T x9 = x*x8;
T x10 = x3*z;
T x11 = x3*x8;
T x12 = -pz[0]*x11 - pz[1]*x9 + pz[2]*x*x8 + pz[3]*x3*x8 - pz[4]*x10 - pz[5]*x7 + pz[6]*x*z + pz[7]*x3*z;
T x13 = x12*x6;
T x14 = -py[0]*x11 - py[1]*x9 + py[2]*x*x8 + py[3]*x3*x8 - py[4]*x10 - py[5]*x7 + py[6]*x*z + py[7]*x3*z;
T x15 = -pz[0]*x5 - pz[1]*x2 - pz[2]*x0 - pz[3]*x4 + pz[4]*x1*x3 + pz[5]*x*x1 + pz[6]*x*y + pz[7]*x3*y;
T x16 = y*z;
T x17 = x8*y;
T x18 = x1*z;
T x19 = x1*x8;
T x20 = -pz[0]*x19 + pz[1]*x1*x8 + pz[2]*x8*y - pz[3]*x17 - pz[4]*x18 + pz[5]*x1*z + pz[6]*y*z - pz[7]*x16;
T x21 = -px[0]*x5 - px[1]*x2 - px[2]*x0 - px[3]*x4 + px[4]*x1*x3 + px[5]*x*x1 + px[6]*x*y + px[7]*x3*y;
T x22 = x14*x21;
T x23 = -py[0]*x19 + py[1]*x1*x8 + py[2]*x8*y - py[3]*x17 - py[4]*x18 + py[5]*x1*z + py[6]*y*z - py[7]*x16;
T x24 = -px[0]*x11 - px[1]*x9 + px[2]*x*x8 + px[3]*x3*x8 - px[4]*x10 - px[5]*x7 + px[6]*x*z + px[7]*x3*z;
T x25 = x15*x24;
T x26 = -px[0]*x19 + px[1]*x1*x8 + px[2]*x8*y - px[3]*x17 - px[4]*x18 + px[5]*x1*z + px[6]*y*z - px[7]*x16;
T x27 = 1.0/(x12*x21*x23 - x13*x26 + x14*x15*x26 - x20*x22 + x20*x24*x6 - x23*x25);
J_inv[0] = x27*(-x13 + x14*x15);
J_inv[1] = x27*(x12*x21 - x25);
J_inv[2] = x27*(-x22 + x24*x6);
J_inv[3] = x27*(-x15*x23 + x20*x6);
J_inv[4] = x27*(x15*x26 - x20*x21);
J_inv[5] = x27*(x21*x23 - x26*x6);
J_inv[6] = x27*(x12*x23 - x14*x20);
J_inv[7] = x27*(-x12*x26 + x20*x24);
J_inv[8] = x27*(x14*x26 - x23*x24);
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
//	- Result: 3*ADD + 3*ASSIGNMENT + 24*MUL
//	- Subexpressions: 12*MUL + 3*SUB
T x0 = y*z;
T x1 = x*x0;
T x2 = 1.0 - z;
T x3 = x2*y;
T x4 = x*x3;
T x5 = 1.0 - y;
T x6 = x5*z;
T x7 = x*x6;
T x8 = 1.0 - x;
T x9 = x0*x8;
T x10 = x2*x5;
T x11 = x*x10;
T x12 = x3*x8;
T x13 = x6*x8;
T x14 = x10*x8;
tx = px[0]*x14 + px[1]*x11 + px[2]*x4 + px[3]*x12 + px[4]*x13 + px[5]*x7 + px[6]*x1 + px[7]*x9;
ty = py[0]*x14 + py[1]*x11 + py[2]*x4 + py[3]*x12 + py[4]*x13 + py[5]*x7 + py[6]*x1 + py[7]*x9;
tz = pz[0]*x14 + pz[1]*x11 + pz[2]*x4 + pz[3]*x12 + pz[4]*x13 + pz[5]*x7 + pz[6]*x1 + pz[7]*x9;
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
				{
	static constexpr int max_iter = 100;
	GeoT tol = 1e-16;
	bool converged = false;
	GeoT residual_norm;
	for(int i = 0; i < max_iter; ++i) {
		//FLOATING POINT OPS!
//	- Result: 13*ADD + 4*ASSIGNMENT + 27*MUL + 3*POW
//	- Subexpressions: 29*ADD + DIV + 171*MUL + 66*SUB
T x0 = x*y;
T x1 = 1.0 - y;
T x2 = x*x1;
T x3 = 1.0 - x;
T x4 = x3*y;
T x5 = x1*x3;
T x6 = -px[0]*x5 - px[1]*x2 - px[2]*x0 - px[3]*x4 + px[4]*x1*x3 + px[5]*x*x1 + px[6]*x*y + px[7]*x3*y;
T x7 = x*z;
T x8 = 1.0 - z;
T x9 = x*x8;
T x10 = x3*z;
T x11 = x3*x8;
T x12 = -py[0]*x11 - py[1]*x9 + py[2]*x*x8 + py[3]*x3*x8 - py[4]*x10 - py[5]*x7 + py[6]*x*z + py[7]*x3*z;
T x13 = x12*x6;
T x14 = -px[0]*x11 - px[1]*x9 + px[2]*x*x8 + px[3]*x3*x8 - px[4]*x10 - px[5]*x7 + px[6]*x*z + px[7]*x3*z;
T x15 = -py[0]*x5 - py[1]*x2 - py[2]*x0 - py[3]*x4 + py[4]*x1*x3 + py[5]*x*x1 + py[6]*x*y + py[7]*x3*y;
T x16 = y*z;
T x17 = pz[7]*x16;
T x18 = x8*y;
T x19 = pz[3]*x18;
T x20 = x1*z;
T x21 = pz[4]*x20;
T x22 = x1*x8;
T x23 = pz[0]*x22;
T x24 = pz[1]*x1*x8 + pz[2]*x8*y + pz[5]*x1*z + pz[6]*y*z - x17 - x19 - x21 - x23;
T x25 = py[7]*x16;
T x26 = py[3]*x18;
T x27 = py[4]*x20;
T x28 = py[0]*x22;
T x29 = py[1]*x1*x8 + py[2]*x8*y + py[5]*x1*z + py[6]*y*z - x25 - x26 - x27 - x28;
T x30 = -pz[0]*x5 - pz[1]*x2 - pz[2]*x0 - pz[3]*x4 + pz[4]*x1*x3 + pz[5]*x*x1 + pz[6]*x*y + pz[7]*x3*y;
T x31 = x14*x30;
T x32 = px[7]*x16;
T x33 = px[3]*x18;
T x34 = px[4]*x20;
T x35 = px[0]*x22;
T x36 = px[1]*x1*x8 + px[2]*x8*y + px[5]*x1*z + px[6]*y*z - x32 - x33 - x34 - x35;
T x37 = -pz[0]*x11 - pz[1]*x9 + pz[2]*x*x8 + pz[3]*x3*x8 - pz[4]*x10 - pz[5]*x7 + pz[6]*x*z + pz[7]*x3*z;
T x38 = x15*x37;
T x39 = 1.0/(x12*x30*x36 - x13*x24 + x14*x15*x24 - x29*x31 + x29*x37*x6 - x36*x38);
T x40 = -pz[1]*x*x22 - pz[2]*x*x18 - pz[5]*x*x20 - pz[6]*x*x16 + tz - x17*x3 - x19*x3 - x21*x3 - x23*x3;
T x41 = x39*x40;
T x42 = -py[1]*x*x22 - py[2]*x*x18 - py[5]*x*x20 - py[6]*x*x16 + ty - x25*x3 - x26*x3 - x27*x3 - x28*x3;
T x43 = x39*x42;
T x44 = -px[1]*x*x22 - px[2]*x*x18 - px[5]*x*x20 - px[6]*x*x16 + tx - x3*x32 - x3*x33 - x3*x34 - x3*x35;
T x45 = x39*x44;
x = tx + x41*(-x13 + x14*x15) + x43*(-x31 + x37*x6) + x45*(x12*x30 - x38);
y = ty + x41*(-x15*x36 + x29*x6) + x43*(-x24*x6 + x30*x36) + x45*(x15*x24 - x29*x30);
z = tz + x41*(x12*x36 - x14*x29) + x43*(x14*x24 - x36*x37) + x45*(-x12*x24 + x29*x37);
residual_norm = pow(x40, 2) + pow(x42, 2) + pow(x44, 2);

		// FIXME not SIMD friendly
		if(residual_norm < tol) {
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
				Result *UTOPIA_RESTRICT gz)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 24*ADD + 24*ASSIGNMENT + 42*MUL
//	- Subexpressions: 29*ADD + DIV + 190*MUL + 5*NEG + 51*SUB
T x0 = x*y;
T x1 = y - 1.0;
T x2 = -x1;
T x3 = x*x2;
T x4 = x - 1.0;
T x5 = -x4;
T x6 = x5*y;
T x7 = x2*x5;
T x8 = -px[0]*x7 - px[1]*x3 - px[2]*x0 - px[3]*x6 + px[4]*x2*x5 + px[5]*x*x2 + px[6]*x*y + px[7]*x5*y;
T x9 = x*z;
T x10 = 1.0 - z;
T x11 = x*x10;
T x12 = x5*z;
T x13 = x10*x5;
T x14 = -pz[0]*x13 - pz[1]*x11 + pz[2]*x*x10 + pz[3]*x10*x5 - pz[4]*x12 - pz[5]*x9 + pz[6]*x*z + pz[7]*x5*z;
T x15 = -px[0]*x13 - px[1]*x11 + px[2]*x*x10 + px[3]*x10*x5 - px[4]*x12 - px[5]*x9 + px[6]*x*z + px[7]*x5*z;
T x16 = -pz[0]*x7 - pz[1]*x3 - pz[2]*x0 - pz[3]*x6 + pz[4]*x2*x5 + pz[5]*x*x2 + pz[6]*x*y + pz[7]*x5*y;
T x17 = x15*x16;
T x18 = x14*x8 - x17;
T x19 = y*z;
T x20 = x10*y;
T x21 = x2*z;
T x22 = x10*x2;
T x23 = -pz[0]*x22 + pz[1]*x10*x2 + pz[2]*x10*y - pz[3]*x20 - pz[4]*x21 + pz[5]*x2*z + pz[6]*y*z - pz[7]*x19;
T x24 = -py[0]*x13 - py[1]*x11 + py[2]*x*x10 + py[3]*x10*x5 - py[4]*x12 - py[5]*x9 + py[6]*x*z + py[7]*x5*z;
T x25 = x24*x8;
T x26 = -py[0]*x22 + py[1]*x10*x2 + py[2]*x10*y - py[3]*x20 - py[4]*x21 + py[5]*x2*z + py[6]*y*z - py[7]*x19;
T x27 = -px[0]*x22 + px[1]*x10*x2 + px[2]*x10*y - px[3]*x20 - px[4]*x21 + px[5]*x2*z + px[6]*y*z - px[7]*x19;
T x28 = -py[0]*x7 - py[1]*x3 - py[2]*x0 - py[3]*x6 + py[4]*x2*x5 + py[5]*x*x2 + py[6]*x*y + py[7]*x5*y;
T x29 = x14*x28;
T x30 = 1.0/(x14*x26*x8 + x15*x23*x28 + x16*x24*x27 - x17*x26 - x23*x25 - x27*x29);
T x31 = x10*x30;
T x32 = x31*x4;
T x33 = x16*x24 - x29;
T x34 = x1*x31;
T x35 = x15*x28 - x25;
T x36 = x30*x7;
T x37 = -x35*x36;
T x38 = x16*x27 - x23*x8;
T x39 = -x16*x26 + x23*x28;
T x40 = x26*x8 - x27*x28;
T x41 = -x36*x40;
T x42 = -x14*x27 + x15*x23;
T x43 = x14*x26 - x23*x24;
T x44 = -x15*x26 + x24*x27;
T x45 = -x36*x44;
T x46 = x3*x30;
T x47 = x35*x46;
T x48 = x11*x30;
T x49 = x18*x48;
T x50 = x40*x46;
T x51 = x38*x48;
T x52 = x44*x46;
T x53 = x42*x48;
T x54 = x0*x30;
T x55 = x35*x54;
T x56 = x20*x30;
T x57 = x33*x56;
T x58 = x40*x54;
T x59 = x39*x56;
T x60 = x44*x54;
T x61 = x43*x56;
T x62 = x30*x6;
T x63 = x35*x62;
T x64 = x40*x62;
T x65 = x44*x62;
T x66 = x12*x30;
T x67 = x18*x66;
T x68 = x21*x30;
T x69 = x33*x68;
T x70 = x38*x66;
T x71 = x39*x68;
T x72 = x42*x66;
T x73 = x43*x68;
T x74 = x30*x9;
T x75 = x18*x74;
T x76 = x38*x74;
T x77 = x42*x74;
T x78 = x19*x30;
T x79 = x33*x78;
T x80 = x39*x78;
T x81 = x43*x78;
gx[0] = x18*x32 + x33*x34 + x37;
gy[0] = x32*x38 + x34*x39 + x41;
gz[0] = x32*x42 + x34*x43 + x45;
gx[1] = x10*x2*x30*x33 - x47 - x49;
gy[1] = x10*x2*x30*x39 - x50 - x51;
gz[1] = x10*x2*x30*x43 - x52 - x53;
gx[2] = x49 - x55 + x57;
gy[2] = x51 - x58 + x59;
gz[2] = x53 - x60 + x61;
gx[3] = x10*x18*x30*x5 - x57 - x63;
gy[3] = x10*x30*x38*x5 - x59 - x64;
gz[3] = x10*x30*x42*x5 - x61 - x65;
gx[4] = -x37 - x67 - x69;
gy[4] = -x41 - x70 - x71;
gz[4] = -x45 - x72 - x73;
gx[5] = x47 + x69 - x75;
gy[5] = x50 + x71 - x76;
gz[5] = x52 + x73 - x77;
gx[6] = x55 + x75 + x79;
gy[6] = x58 + x76 + x80;
gz[6] = x60 + x77 + x81;
gx[7] = x63 + x67 - x79;
gy[7] = x64 + x70 - x80;
gz[7] = x65 + x72 - x81;
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
//	- Result: 24*ADD + 33*ASSIGNMENT + 50*MUL
//	- Subexpressions: 29*ADD + DIV + 190*MUL + 5*NEG + 51*SUB
T x0 = x - 1.0;
T x1 = -x0;
T x2 = y - 1.0;
T x3 = -x2;
T x4 = 1.0 - z;
T x5 = x3*x4;
T x6 = x4*y;
T x7 = x3*z;
T x8 = y*z;
T x9 = -pz[0]*x5 + pz[1]*x3*x4 + pz[2]*x4*y - pz[3]*x6 - pz[4]*x7 + pz[5]*x3*z + pz[6]*y*z - pz[7]*x8;
T x10 = x*y;
T x11 = x*x3;
T x12 = x1*y;
T x13 = x1*x3;
T x14 = -px[0]*x13 - px[1]*x11 - px[2]*x10 - px[3]*x12 + px[4]*x1*x3 + px[5]*x*x3 + px[6]*x*y + px[7]*x1*y;
T x15 = x*z;
T x16 = x*x4;
T x17 = x1*z;
T x18 = x1*x4;
T x19 = -py[0]*x18 - py[1]*x16 + py[2]*x*x4 + py[3]*x1*x4 - py[4]*x17 - py[5]*x15 + py[6]*x*z + py[7]*x1*z;
T x20 = x14*x19;
T x21 = -py[0]*x5 + py[1]*x3*x4 + py[2]*x4*y - py[3]*x6 - py[4]*x7 + py[5]*x3*z + py[6]*y*z - py[7]*x8;
T x22 = -px[0]*x18 - px[1]*x16 + px[2]*x*x4 + px[3]*x1*x4 - px[4]*x17 - px[5]*x15 + px[6]*x*z + px[7]*x1*z;
T x23 = -pz[0]*x13 - pz[1]*x11 - pz[2]*x10 - pz[3]*x12 + pz[4]*x1*x3 + pz[5]*x*x3 + pz[6]*x*y + pz[7]*x1*y;
T x24 = x22*x23;
T x25 = -px[0]*x5 + px[1]*x3*x4 + px[2]*x4*y - px[3]*x6 - px[4]*x7 + px[5]*x3*z + px[6]*y*z - px[7]*x8;
T x26 = -py[0]*x13 - py[1]*x11 - py[2]*x10 - py[3]*x12 + py[4]*x1*x3 + py[5]*x*x3 + py[6]*x*y + py[7]*x1*y;
T x27 = -pz[0]*x18 - pz[1]*x16 + pz[2]*x*x4 + pz[3]*x1*x4 - pz[4]*x17 - pz[5]*x15 + pz[6]*x*z + pz[7]*x1*z;
T x28 = x26*x27;
T x29 = x14*x21*x27 + x19*x23*x25 - x20*x9 - x21*x24 + x22*x26*x9 - x25*x28;
T x30 = x14*x27 - x24;
T x31 = 1.0/x29;
T x32 = x31*x4;
T x33 = x0*x32;
T x34 = x19*x23 - x28;
T x35 = x2*x32;
T x36 = -x20 + x22*x26;
T x37 = x13*x31;
T x38 = -x36*x37;
T x39 = -x14*x9 + x23*x25;
T x40 = -x21*x23 + x26*x9;
T x41 = x14*x21 - x25*x26;
T x42 = -x37*x41;
T x43 = x22*x9 - x25*x27;
T x44 = -x19*x9 + x21*x27;
T x45 = x19*x25 - x21*x22;
T x46 = -x37*x45;
T x47 = x11*x31;
T x48 = x36*x47;
T x49 = x16*x31;
T x50 = x30*x49;
T x51 = x41*x47;
T x52 = x39*x49;
T x53 = x45*x47;
T x54 = x43*x49;
T x55 = x10*x31;
T x56 = x36*x55;
T x57 = x31*x6;
T x58 = x34*x57;
T x59 = x41*x55;
T x60 = x40*x57;
T x61 = x45*x55;
T x62 = x44*x57;
T x63 = x12*x31;
T x64 = x36*x63;
T x65 = x41*x63;
T x66 = x45*x63;
T x67 = x17*x31;
T x68 = x30*x67;
T x69 = x31*x7;
T x70 = x34*x69;
T x71 = x39*x67;
T x72 = x40*x69;
T x73 = x43*x67;
T x74 = x44*x69;
T x75 = x15*x31;
T x76 = x30*x75;
T x77 = x39*x75;
T x78 = x43*x75;
T x79 = x31*x8;
T x80 = x34*x79;
T x81 = x40*x79;
T x82 = x44*x79;
f[0] = x1*x5;
f[1] = x*x5;
f[2] = x*x6;
f[3] = x1*x6;
f[4] = x1*x7;
f[5] = x*x7;
f[6] = x*x8;
f[7] = x1*x8;
measure_value = x29;
gx[0] = x30*x33 + x34*x35 + x38;
gy[0] = x33*x39 + x35*x40 + x42;
gz[0] = x33*x43 + x35*x44 + x46;
gx[1] = x3*x31*x34*x4 - x48 - x50;
gy[1] = x3*x31*x4*x40 - x51 - x52;
gz[1] = x3*x31*x4*x44 - x53 - x54;
gx[2] = x50 - x56 + x58;
gy[2] = x52 - x59 + x60;
gz[2] = x54 - x61 + x62;
gx[3] = x1*x30*x31*x4 - x58 - x64;
gy[3] = x1*x31*x39*x4 - x60 - x65;
gz[3] = x1*x31*x4*x43 - x62 - x66;
gx[4] = -x38 - x68 - x70;
gy[4] = -x42 - x71 - x72;
gz[4] = -x46 - x73 - x74;
gx[5] = x48 + x70 - x76;
gy[5] = x51 + x72 - x77;
gz[5] = x53 + x74 - x78;
gx[6] = x56 + x76 + x80;
gy[6] = x59 + x77 + x81;
gz[6] = x61 + x78 + x82;
gx[7] = x64 + x68 - x80;
gy[7] = x65 + x71 - x81;
gz[7] = x66 + x73 - x82;
			}


		};
	}
}

#endif // UTOPIA_TPL_FE_Hex8_3_IMPL_hpp
