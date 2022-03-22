#ifndef UTOPIA_TPL_FE_Pentatope5_4_IMPL_hpp
#define UTOPIA_TPL_FE_Pentatope5_4_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"

// #include "utopia_fe_Pentatope5.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of Pentatope5 for dimension 4
		 */
		template<typename T, typename GeoT>
		class Pentatope5 {
		public:
			static constexpr int Dim = 4;
			static constexpr int NNodes = 5;
			static constexpr int Order = 1;

			using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

			UTOPIA_FUNCTION static constexpr const char* class_name() { return "Pentatope5"; }

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
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t)
			{
				T measure_value;
				
// Unused variables
UTOPIA_UNUSED(x);
UTOPIA_UNUSED(y);
UTOPIA_UNUSED(z);
UTOPIA_UNUSED(t);
//FLOATING POINT OPS!
//	- Result: 9*ADD + ASSIGNMENT + 32*MUL
//	- Subexpressions: 11*MUL + NEG + 12*SUB
T x0 = px[0] - px[2];
T x1 = -x0;
T x2 = -py[0] + py[4];
T x3 = -pz[0] + pz[3];
T x4 = x2*x3;
T x5 = -px[0] + px[3];
T x6 = -pz[0] + pz[4];
T x7 = -py[0] + py[2];
T x8 = x6*x7;
T x9 = -px[0] + px[4];
T x10 = -py[0] + py[3];
T x11 = -pz[0] + pz[2];
T x12 = x10*x11;
T x13 = x10*x6;
T x14 = x11*x2;
T x15 = x3*x7;
T x16 = -pt[0] + pt[2];
T x17 = x16*x9;
T x18 = -pt[0] + pt[3];
T x19 = x1*x18;
T x20 = -pt[0] + pt[4];
T x21 = x20*x5;
T x22 = x16*x5;
T x23 = x18*x9;
measure_value = (-pt[0] + pt[1])*(x0*x13 + x1*x4 + x12*x9 - x14*x5 - x15*x9 + x5*x8) + (-px[0] + px[1])*(-x12*x20 + x13*x16 + x14*x18 + x15*x20 - x16*x4 - x18*x8) + (-py[0] + py[1])*(x0*x20*x3 + x11*x21 - x11*x23 + x17*x3 + x19*x6 - x22*x6) + (-pz[0] + pz[1])*(x1*x10*x20 - x10*x17 - x19*x2 + x2*x22 - x21*x7 + x23*x7);
				return measure_value;
			}

			UTOPIA_FUNCTION static void jacobian(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				GeoT *UTOPIA_RESTRICT J)
			{
				using namespace utopia::device;
				// Automatically generated
				
// Unused variables
UTOPIA_UNUSED(x);
UTOPIA_UNUSED(y);
UTOPIA_UNUSED(z);
UTOPIA_UNUSED(t);
//FLOATING POINT OPS!
//	- Result: 16*ADD + 16*ASSIGNMENT + 16*MUL
//	- Subexpressions: 0
J[0] = -px[0] + px[1];
J[1] = -py[0] + py[1];
J[2] = -pz[0] + pz[1];
J[3] = -pt[0] + pt[1];
J[4] = -px[0] + px[2];
J[5] = -py[0] + py[2];
J[6] = -pz[0] + pz[2];
J[7] = -pt[0] + pt[2];
J[8] = -px[0] + px[3];
J[9] = -py[0] + py[3];
J[10] = -pz[0] + pz[3];
J[11] = -pt[0] + pt[3];
J[12] = -px[0] + px[4];
J[13] = -py[0] + py[4];
J[14] = -pz[0] + pz[4];
J[15] = -pt[0] + pt[4];
			}

			UTOPIA_FUNCTION static void jacobian_inverse(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				GeoT *UTOPIA_RESTRICT J_inv)
			{
				using namespace utopia::device;
				// Automatically generated
				
// Unused variables
UTOPIA_UNUSED(x);
UTOPIA_UNUSED(y);
UTOPIA_UNUSED(z);
UTOPIA_UNUSED(t);
//FLOATING POINT OPS!
//	- Result: 11*ADD + 16*ASSIGNMENT + 81*MUL
//	- Subexpressions: 13*ADD + DIV + 56*MUL + 3*NEG + 26*SUB
T x0 = -pt[0] + pt[2];
T x1 = -py[0] + py[3];
T x2 = -pz[0] + pz[4];
T x3 = x1*x2;
T x4 = -pt[0] + pt[3];
T x5 = -py[0] + py[4];
T x6 = -pz[0] + pz[2];
T x7 = x5*x6;
T x8 = -pt[0] + pt[4];
T x9 = -pz[0] + pz[3];
T x10 = -py[0] + py[2];
T x11 = x10*x9;
T x12 = x5*x9;
T x13 = x10*x2;
T x14 = x1*x6;
T x15 = -x0*x12 + x0*x3 + x11*x8 - x13*x4 - x14*x8 + x4*x7;
T x16 = -pt[0] + pt[1];
T x17 = px[0] - px[2];
T x18 = -x17;
T x19 = -px[0] + px[3];
T x20 = -px[0] + px[4];
T x21 = -x11*x20 + x12*x18 + x13*x19 + x14*x20 + x17*x3 - x19*x7;
T x22 = py[0] - py[1];
T x23 = -x22;
T x24 = x0*x20;
T x25 = x18*x4;
T x26 = x19*x8;
T x27 = x8*x9;
T x28 = x0*x19;
T x29 = x20*x4;
T x30 = x17*x27 + x2*x25 - x2*x28 + x24*x9 + x26*x6 - x29*x6;
T x31 = px[0] - px[1];
T x32 = -x31;
T x33 = -pz[0] + pz[1];
T x34 = x1*x8;
T x35 = -x1*x24 - x10*x26 + x10*x29 + x18*x34 - x25*x5 + x28*x5;
T x36 = 1.0/(x15*x32 + x16*x21 + x23*x30 + x33*x35);
T x37 = x2*x4;
T x38 = x33*x4;
T x39 = x0*x9;
T x40 = x4*x6;
T x41 = x0*x1;
T x42 = x16*x2;
T x43 = x16*x20;
T x44 = x0*x32;
T x45 = x18*x8;
T x46 = x16*x18;
T x47 = x16*x19;
T x48 = x33*x5;
T x49 = x20*x23;
T x50 = x2*x23;
T x51 = x20*x33;
J_inv[0] = x15*x36;
J_inv[1] = x36*(x12*x16 - x16*x3 + x22*x27 + x23*x37 + x33*x34 - x38*x5);
J_inv[2] = 0;
J_inv[3] = x36*(x10*x38 - x11*x16 + x14*x16 + x22*x40 + x23*x39 - x33*x41);
J_inv[4] = x30*x36;
J_inv[5] = x36*(x19*x42 - x26*x33 + x27*x32 + x29*x33 - x32*x37 - x43*x9);
J_inv[6] = x36*(-x18*x42 + x2*x44 - x24*x33 + x31*x6*x8 + x33*x45 + x43*x6);
J_inv[7] = x36*(-x25*x33 + x28*x33 - x32*x39 + x32*x40 + x46*x9 - x47*x6);
J_inv[8] = x35*x36;
J_inv[9] = x36*(x1*x43 + x23*x26 - x23*x29 + x31*x34 + x32*x4*x5 - x47*x5);
J_inv[10] = x36*(x10*x32*x8 - x10*x43 + x23*x24 - x23*x45 - x44*x5 + x46*x5);
J_inv[11] = x36*(-x1*x46 + x10*x31*x4 + x10*x47 + x23*x25 - x23*x28 + x32*x41);
J_inv[12] = x21*x36;
J_inv[13] = x36*(-x1*x51 - x12*x32 + x19*x48 - x19*x50 + x3*x32 + x49*x9);
J_inv[14] = x36*(x10*x51 + x13*x31 - x18*x48 + x18*x50 + x32*x7 - x49*x6);
J_inv[15] = x36*(x1*x18*x33 - x10*x19*x33 + x11*x32 - x14*x32 - x18*x23*x9 + x19*x23*x6);
			}

			UTOPIA_FUNCTION static void transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				GeoT &tx,
				GeoT &ty,
				GeoT &tz,
				GeoT &tt)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 4*ADD + 4*ASSIGNMENT + 20*MUL
//	- Subexpressions: 4*SUB
T x0 = -t - x - y - z + 1;
tx = px[0]*x0 + px[1]*x + px[2]*y + px[3]*z + px[4]*t;
ty = py[0]*x0 + py[1]*x + py[2]*y + py[3]*z + py[4]*t;
tz = pz[0]*x0 + pz[1]*x + pz[2]*y + pz[3]*z + pz[4]*t;
tt = pt[0]*x0 + pt[1]*x + pt[2]*y + pt[3]*z + pt[4]*t;
			}

			UTOPIA_FUNCTION static void inverse_transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T tx,
				const T ty,
				const T tz,
				const T tt,
				GeoT &x,
				GeoT &y,
				GeoT &z,
				GeoT &t)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 4*ADD + 4*ASSIGNMENT + 200*MUL
//	- Subexpressions: 103*ADD + DIV + 324*MUL + 108*SUB
T x0 = py[4]*pz[2];
T x1 = px[0]*tt;
T x2 = x0*x1;
T x3 = py[0]*pz[4];
T x4 = px[2]*tt;
T x5 = x3*x4;
T x6 = py[4]*pz[3];
T x7 = py[2]*pz[4];
T x8 = px[3]*tt;
T x9 = py[2]*pz[0];
T x10 = px[4]*tt;
T x11 = x10*x9;
T x12 = py[3]*pz[2];
T x13 = pt[0]*tx;
T x14 = x13*x7;
T x15 = py[3]*pz[4];
T x16 = pt[2]*tx;
T x17 = py[4]*pz[0];
T x18 = x16*x17;
T x19 = pt[3]*tx;
T x20 = py[0]*pz[2];
T x21 = pt[4]*tx;
T x22 = x20*x21;
T x23 = py[2]*pz[3];
T x24 = pt[0]*px[4];
T x25 = pz[2]*ty;
T x26 = x24*x25;
T x27 = pt[2]*px[0];
T x28 = pz[4]*ty;
T x29 = x27*x28;
T x30 = pt[2]*px[4];
T x31 = pz[3]*ty;
T x32 = pt[3]*px[2];
T x33 = pt[4]*px[2];
T x34 = pz[0]*ty;
T x35 = x33*x34;
T x36 = pt[4]*px[3];
T x37 = pt[0]*px[2];
T x38 = py[4]*tz;
T x39 = x37*x38;
T x40 = pt[2]*px[3];
T x41 = py[0]*tz;
T x42 = x30*x41;
T x43 = pt[3]*px[4];
T x44 = py[2]*tz;
T x45 = pt[4]*px[0];
T x46 = x44*x45;
T x47 = py[3]*tz;
T x48 = x15*x37;
T x49 = pt[0]*px[3];
T x50 = x0*x49;
T x51 = x23*x24;
T x52 = x27*x6;
T x53 = x3*x40;
T x54 = py[3]*pz[0];
T x55 = x30*x54;
T x56 = pt[3]*px[0];
T x57 = x56*x7;
T x58 = x17*x32;
T x59 = x20*x43;
T x60 = x12*x45;
T x61 = py[0]*pz[3];
T x62 = x33*x61;
T x63 = x36*x9;
T x64 = x1*x7;
T x65 = x17*x4;
T x66 = x10*x20;
T x67 = x0*x13;
T x68 = x16*x3;
T x69 = x21*x9;
T x70 = x28*x37;
T x71 = x30*x34;
T x72 = x25*x45;
T x73 = x24*x44;
T x74 = x27*x38;
T x75 = x33*x41;
T x76 = x1*x15 - x1*x6 - x10*x54 + x10*x61 - x13*x15 + x13*x6 - x17*x19 + x17*x8 + x19*x3 + x21*x54 - x21*x61 - x24*x31 + x24*x47 + x28*x49 - x28*x56 - x3*x8 + x31*x45 - x34*x36 + x34*x43 + x36*x41 - x38*x49 + x38*x56 - x41*x43 - x45*x47;
T x77 = -x1*x12 + x1*x23 + x12*x13 - x13*x23 - x16*x54 + x16*x61 - x19*x20 + x19*x9 + x20*x8 - x25*x49 + x25*x56 - x27*x31 + x27*x47 + x31*x37 - x32*x34 + x32*x41 + x34*x40 - x37*x47 + x4*x54 - x4*x61 - x40*x41 + x44*x49 - x44*x56 - x8*x9;
T x78 = pt[0]*px[1];
T x79 = py[1]*pz[4];
T x80 = py[2]*pz[1];
T x81 = pt[1]*px[0];
T x82 = pt[1]*px[2];
T x83 = pt[1]*px[3];
T x84 = pt[1]*px[4];
T x85 = py[4]*pz[1];
T x86 = pt[2]*px[1];
T x87 = py[1]*pz[0];
T x88 = py[3]*pz[1];
T x89 = pt[3]*px[1];
T x90 = py[1]*pz[2];
T x91 = pt[4]*px[1];
T x92 = py[0]*pz[1];
T x93 = py[1]*pz[3];
T x94 = x7*x78;
T x95 = x37*x85;
T x96 = x24*x90;
T x97 = x0*x81;
T x98 = x3*x82;
T x99 = x84*x9;
T x100 = x27*x79;
T x101 = x17*x86;
T x102 = x30*x92;
T x103 = x45*x80;
T x104 = x20*x91;
T x105 = x33*x87;
T x106 = x15*x78 - x15*x81 - x17*x83 + x17*x89 - x24*x88 + x24*x93 + x3*x83 - x3*x89 + x36*x87 - x36*x92 - x43*x87 + x43*x92 + x45*x88 - x45*x93 - x49*x79 + x49*x85 + x54*x84 - x54*x91 + x56*x79 - x56*x85 - x6*x78 + x6*x81 - x61*x84 + x61*x91;
T x107 = -x12*x78 + x12*x81 - x20*x83 + x20*x89 + x23*x78 - x23*x81 - x27*x88 + x27*x93 + x32*x87 - x32*x92 + x37*x88 - x37*x93 - x40*x87 + x40*x92 - x49*x80 + x49*x90 - x54*x82 + x54*x86 + x56*x80 - x56*x90 + x61*x82 - x61*x86 + x83*x9 - x89*x9;
T x108 = 1.0/(x0*x56 + x0*x78 + x0*x83 - x0*x89 - x100 - x101 - x102 - x103 - x104 - x105 + x106 + x107 + x12*x24 - x12*x84 + x12*x91 + x15*x27 + x15*x82 - x15*x86 + x17*x40 + x17*x82 + x20*x36 + x20*x84 + x23*x45 + x23*x84 - x23*x91 + x24*x80 + x27*x85 + x3*x32 + x3*x86 + x30*x61 + x30*x87 + x30*x88 - x30*x93 - x32*x79 + x32*x85 + x33*x54 - x33*x88 + x33*x92 + x33*x93 + x36*x80 - x36*x90 + x37*x6 + x37*x79 + x40*x79 - x40*x85 - x43*x80 + x43*x9 + x43*x90 + x45*x90 - x48 + x49*x7 - x50 - x51 - x52 - x53 - x55 - x57 - x58 - x59 - x6*x82 + x6*x86 - x60 - x62 - x63 + x7*x81 - x7*x83 + x7*x89 + x9*x91 - x94 - x95 - x96 - x97 - x98 - x99);
T x109 = x1*x93;
T x110 = px[1]*tt;
T x111 = x110*x54;
T x112 = x8*x92;
T x113 = x13*x88;
T x114 = pt[1]*tx;
T x115 = x114*x61;
T x116 = x19*x87;
T x117 = x31*x78;
T x118 = x34*x83;
T x119 = pz[1]*ty;
T x120 = x119*x56;
T x121 = py[1]*tz;
T x122 = x121*x49;
T x123 = x47*x81;
T x124 = x41*x89;
T x125 = x1*x88;
T x126 = x110*x61;
T x127 = x8*x87;
T x128 = x13*x93;
T x129 = x114*x54;
T x130 = x19*x92;
T x131 = x119*x49;
T x132 = x31*x81;
T x133 = x34*x89;
T x134 = x47*x78;
T x135 = x41*x83;
T x136 = x121*x56;
T x137 = -x1*x79 + x1*x85 + x10*x87 - x10*x92 - x110*x17 + x110*x3 + x114*x17 - x114*x3 + x119*x24 - x119*x45 - x121*x24 + x121*x45 + x13*x79 - x13*x85 - x21*x87 + x21*x92 - x28*x78 + x28*x81 - x34*x84 + x34*x91 + x38*x78 - x38*x81 + x41*x84 - x41*x91;
T x138 = -x1*x80 + x1*x90 - x110*x20 + x110*x9 + x114*x20 - x114*x9 + x119*x27 - x119*x37 - x121*x27 + x121*x37 + x13*x80 - x13*x90 + x16*x87 - x16*x92 + x25*x78 - x25*x81 + x34*x82 - x34*x86 - x4*x87 + x4*x92 - x41*x82 + x41*x86 - x44*x78 + x44*x81;
x = x108*(pt[0]*px[2]*py[4]*pz[3] + pt[0]*px[3]*py[2]*pz[4] + pt[0]*px[4]*py[3]*pz[2] + pt[2]*px[0]*py[3]*pz[4] + pt[2]*px[3]*py[4]*pz[0] + pt[2]*px[3]*pz[4]*ty + pt[2]*px[4]*py[0]*pz[3] + pt[2]*px[4]*py[3]*tz + pt[2]*py[4]*pz[3]*tx + pt[3]*px[0]*py[4]*pz[2] + pt[3]*px[2]*py[0]*pz[4] + pt[3]*px[2]*py[4]*tz + pt[3]*px[4]*py[2]*pz[0] + pt[3]*px[4]*pz[2]*ty + pt[3]*py[2]*pz[4]*tx + pt[4]*px[0]*py[2]*pz[3] + pt[4]*px[2]*py[3]*pz[0] + pt[4]*px[2]*pz[3]*ty + pt[4]*px[3]*py[0]*pz[2] + pt[4]*px[3]*py[2]*tz + pt[4]*py[3]*pz[2]*tx + px[2]*py[3]*pz[4]*tt + px[3]*py[4]*pz[2]*tt + px[4]*py[2]*pz[3]*tt - x0*x19 - x10*x12 - x11 - x14 - x15*x16 - x18 - x2 - x21*x23 - x22 - x25*x36 - x26 - x28*x32 - x29 - x30*x31 - x33*x47 - x35 - x38*x40 - x39 - x4*x6 - x42 - x43*x44 - x46 - x48 - x5 - x50 - x51 - x52 - x53 - x55 - x57 - x58 - x59 - x60 - x62 - x63 + x64 + x65 + x66 + x67 + x68 + x69 - x7*x8 + x70 + x71 + x72 + x73 + x74 + x75 - x76 - x77);
y = x108*(x10*x88 - x10*x93 + x106 + x109 - x110*x15 + x110*x6 + x111 + x112 + x113 + x114*x15 - x114*x6 + x115 + x116 + x117 + x118 + x119*x36 - x119*x43 + x120 - x121*x36 + x121*x43 + x122 + x123 + x124 - x125 - x126 - x127 - x128 - x129 - x130 - x131 - x132 - x133 - x134 - x135 - x136 + x137 - x19*x79 + x19*x85 - x21*x88 + x21*x93 - x28*x83 + x28*x89 + x31*x84 - x31*x91 + x38*x83 - x38*x89 - x47*x84 + x47*x91 + x76 + x79*x8 - x8*x85);
z = x108*(pt[0]*px[1]*py[4]*pz[2] + pt[0]*px[2]*py[1]*pz[4] + pt[0]*px[4]*py[2]*pz[1] + pt[1]*px[0]*py[2]*pz[4] + pt[1]*px[2]*py[4]*pz[0] + pt[1]*px[2]*pz[4]*ty + pt[1]*px[4]*py[0]*pz[2] + pt[1]*px[4]*py[2]*tz + pt[1]*py[4]*pz[2]*tx + pt[2]*px[0]*py[4]*pz[1] + pt[2]*px[1]*py[0]*pz[4] + pt[2]*px[1]*py[4]*tz + pt[2]*px[4]*py[1]*pz[0] + pt[2]*px[4]*pz[1]*ty + pt[2]*py[1]*pz[4]*tx + pt[4]*px[0]*py[1]*pz[2] + pt[4]*px[1]*py[2]*pz[0] + pt[4]*px[1]*pz[2]*ty + pt[4]*px[2]*py[0]*pz[1] + pt[4]*px[2]*py[1]*tz + pt[4]*py[2]*pz[1]*tx + px[1]*py[2]*pz[4]*tt + px[2]*py[4]*pz[1]*tt + px[4]*py[1]*pz[2]*tt - x0*x110 - x10*x80 - x100 - x101 - x102 - x103 - x104 - x105 + x11 - x114*x7 - x119*x33 - x121*x30 - x137 - x138 + x14 - x16*x85 + x18 + x2 - x21*x90 + x22 - x25*x84 + x26 - x28*x86 + x29 + x35 - x38*x82 + x39 - x4*x79 + x42 - x44*x91 + x46 + x5 - x64 - x65 - x66 - x67 - x68 - x69 - x70 - x71 - x72 - x73 - x74 - x75 - x94 - x95 - x96 - x97 - x98 - x99);
t = x108*(x107 - x109 + x110*x12 - x110*x23 - x111 - x112 - x113 - x114*x12 + x114*x23 - x115 - x116 - x117 - x118 + x119*x32 - x119*x40 - x120 - x121*x32 + x121*x40 - x122 - x123 - x124 + x125 + x126 + x127 + x128 + x129 + x130 + x131 + x132 + x133 + x134 + x135 + x136 + x138 + x16*x88 - x16*x93 - x19*x80 + x19*x90 + x25*x83 - x25*x89 - x31*x82 + x31*x86 - x4*x88 + x4*x93 - x44*x83 + x44*x89 + x47*x82 - x47*x86 + x77 + x8*x80 - x8*x90);
			}

			UTOPIA_FUNCTION static void gradient(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				// Output
				Result *UTOPIA_RESTRICT gx,
				Result *UTOPIA_RESTRICT gy,
				Result *UTOPIA_RESTRICT gz,
				Result *UTOPIA_RESTRICT gt)
			{
				using namespace utopia::device;
				// Automatically generated
				
// Unused variables
UTOPIA_UNUSED(x);
UTOPIA_UNUSED(y);
UTOPIA_UNUSED(z);
UTOPIA_UNUSED(t);
//FLOATING POINT OPS!
//	- Result: 4*ADD + 20*ASSIGNMENT + 15*MUL
//	- Subexpressions: 41*ADD + DIV + 145*MUL + 3*NEG + 53*SUB
T x0 = -pt[0] + pt[1];
T x1 = -py[0] + py[3];
T x2 = -pz[0] + pz[2];
T x3 = x1*x2;
T x4 = py[0] - py[1];
T x5 = -x4;
T x6 = -pz[0] + pz[3];
T x7 = -pt[0] + pt[2];
T x8 = x6*x7;
T x9 = -py[0] + py[2];
T x10 = -pz[0] + pz[1];
T x11 = -pt[0] + pt[3];
T x12 = x10*x11;
T x13 = x11*x2;
T x14 = x6*x9;
T x15 = x1*x7;
T x16 = px[0] - px[2];
T x17 = -x16;
T x18 = -py[0] + py[4];
T x19 = x18*x6;
T x20 = -px[0] + px[3];
T x21 = -pz[0] + pz[4];
T x22 = x21*x9;
T x23 = -px[0] + px[4];
T x24 = x1*x21;
T x25 = x18*x2;
T x26 = -x14*x23 + x16*x24 + x17*x19 + x20*x22 - x20*x25 + x23*x3;
T x27 = x23*x7;
T x28 = x11*x17;
T x29 = -pt[0] + pt[4];
T x30 = x20*x29;
T x31 = x29*x6;
T x32 = x20*x7;
T x33 = x11*x23;
T x34 = x16*x31 + x2*x30 - x2*x33 + x21*x28 - x21*x32 + x27*x6;
T x35 = px[0] - px[1];
T x36 = -x35;
T x37 = -x11*x22 + x11*x25 + x14*x29 - x19*x7 + x24*x7 - x29*x3;
T x38 = x1*x29;
T x39 = -x1*x27 + x17*x38 - x18*x28 + x18*x32 - x30*x9 + x33*x9;
T x40 = 1.0/(x0*x26 + x10*x39 + x34*x5 + x36*x37);
T x41 = x40*(-x0*x14 + x0*x3 - x10*x15 + x12*x9 + x13*x4 + x5*x8);
T x42 = x11*x21;
T x43 = x40*(x0*x19 - x0*x24 + x10*x38 - x12*x18 + x31*x4 + x42*x5);
T x44 = x37*x40;
T x45 = x0*x23;
T x46 = x36*x7;
T x47 = x17*x29;
T x48 = x0*x17;
T x49 = x40*(-x10*x27 + x10*x47 + x2*x29*x35 + x2*x45 + x21*x46 - x21*x48);
T x50 = x34*x40;
T x51 = x0*x6;
T x52 = x0*x20;
T x53 = x40*(-x10*x28 + x10*x32 + x13*x36 + x17*x51 - x2*x52 - x36*x8);
T x54 = x40*(-x10*x30 + x10*x33 + x21*x52 - x23*x51 + x31*x36 - x36*x42);
T x55 = x40*(-x1*x48 + x11*x35*x9 + x15*x36 + x28*x5 - x32*x5 + x52*x9);
T x56 = x40*(x1*x45 + x11*x18*x36 - x18*x52 + x30*x5 - x33*x5 + x35*x38);
T x57 = x40*(-x18*x46 + x18*x48 + x27*x5 + x29*x36*x9 - x45*x9 - x47*x5);
T x58 = x39*x40;
T x59 = x17*x5;
T x60 = x10*x9;
T x61 = x10*x17;
T x62 = x2*x5;
T x63 = x40*(-x18*x61 + x21*x59 + x22*x35 + x23*x60 - x23*x62 + x25*x36);
T x64 = x26*x40;
T x65 = x40*(x1*x61 + x14*x36 - x20*x60 + x20*x62 - x3*x36 - x59*x6);
T x66 = x40*(-x1*x10*x23 + x10*x18*x20 - x19*x36 - x20*x21*x5 + x23*x5*x6 + x24*x36);
gx[0] = -x41 - x43 - x44;
gy[0] = -x49 - x50 - x53 - x54;
gz[0] = -x55 - x56 - x57 - x58;
gt[0] = -x63 - x64 - x65 - x66;
gx[1] = x44;
gy[1] = x50;
gz[1] = x58;
gt[1] = x64;
gx[2] = x43;
gy[2] = x54;
gz[2] = x56;
gt[2] = x66;
gx[3] = 0;
gy[3] = x49;
gz[3] = x57;
gt[3] = x63;
gx[4] = x41;
gy[4] = x53;
gz[4] = x55;
gt[4] = x65;
			}

			UTOPIA_FUNCTION static void value(
				const T x,
				const T y,
				const T z,
				const T t,
				Result *UTOPIA_RESTRICT f
				)
			{
				using namespace utopia::device;
			    // Automatically generated
				//FLOATING POINT OPS!
//	- Result: ADD + 5*ASSIGNMENT + 4*MUL
//	- Subexpressions: 0
f[0] = -t - x - y - z + 1;
f[1] = x;
f[2] = y;
f[3] = z;
f[4] = t;
			}

			UTOPIA_FUNCTION static void eval(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				// Output
				Result *UTOPIA_RESTRICT f,
				Result *UTOPIA_RESTRICT gx,
				Result *UTOPIA_RESTRICT gy,
				Result *UTOPIA_RESTRICT gz,
				Result *UTOPIA_RESTRICT gt,
				T &measure_value)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 5*ADD + 26*ASSIGNMENT + 19*MUL
//	- Subexpressions: 41*ADD + DIV + 145*MUL + 3*NEG + 53*SUB
T x0 = -pt[0] + pt[1];
T x1 = px[0] - px[2];
T x2 = -x1;
T x3 = -py[0] + py[4];
T x4 = -pz[0] + pz[3];
T x5 = x3*x4;
T x6 = -px[0] + px[3];
T x7 = -pz[0] + pz[4];
T x8 = -py[0] + py[2];
T x9 = x7*x8;
T x10 = -px[0] + px[4];
T x11 = -py[0] + py[3];
T x12 = -pz[0] + pz[2];
T x13 = x11*x12;
T x14 = x11*x7;
T x15 = x12*x3;
T x16 = x4*x8;
T x17 = x1*x14 + x10*x13 - x10*x16 - x15*x6 + x2*x5 + x6*x9;
T x18 = py[0] - py[1];
T x19 = -x18;
T x20 = -pt[0] + pt[2];
T x21 = x10*x20;
T x22 = -pt[0] + pt[3];
T x23 = x2*x22;
T x24 = -pt[0] + pt[4];
T x25 = x24*x6;
T x26 = x24*x4;
T x27 = x20*x6;
T x28 = x10*x22;
T x29 = x1*x26 + x12*x25 - x12*x28 + x21*x4 + x23*x7 - x27*x7;
T x30 = px[0] - px[1];
T x31 = -x30;
T x32 = -x13*x24 + x14*x20 + x15*x22 + x16*x24 - x20*x5 - x22*x9;
T x33 = -pz[0] + pz[1];
T x34 = x11*x24;
T x35 = -x11*x21 + x2*x34 - x23*x3 - x25*x8 + x27*x3 + x28*x8;
T x36 = x0*x17 + x19*x29 + x31*x32 + x33*x35;
T x37 = x20*x4;
T x38 = x22*x33;
T x39 = x12*x22;
T x40 = x11*x20;
T x41 = 1.0/x36;
T x42 = x41*(x0*x13 - x0*x16 + x18*x39 + x19*x37 - x33*x40 + x38*x8);
T x43 = x22*x7;
T x44 = x41*(-x0*x14 + x0*x5 + x18*x26 + x19*x43 - x3*x38 + x33*x34);
T x45 = x32*x41;
T x46 = x0*x10;
T x47 = x20*x31;
T x48 = x2*x24;
T x49 = x0*x2;
T x50 = x41*(x12*x24*x30 + x12*x46 - x21*x33 + x33*x48 + x47*x7 - x49*x7);
T x51 = x29*x41;
T x52 = x0*x4;
T x53 = x0*x6;
T x54 = x41*(-x12*x53 + x2*x52 - x23*x33 + x27*x33 - x31*x37 + x31*x39);
T x55 = x41*(-x10*x52 - x25*x33 + x26*x31 + x28*x33 - x31*x43 + x53*x7);
T x56 = x41*(-x11*x49 + x19*x23 - x19*x27 + x22*x30*x8 + x31*x40 + x53*x8);
T x57 = x41*(x11*x46 + x19*x25 - x19*x28 + x22*x3*x31 - x3*x53 + x30*x34);
T x58 = x41*(x19*x21 - x19*x48 + x24*x31*x8 - x3*x47 + x3*x49 - x46*x8);
T x59 = x35*x41;
T x60 = x19*x2;
T x61 = x33*x8;
T x62 = x2*x33;
T x63 = x12*x19;
T x64 = x41*(x10*x61 - x10*x63 + x15*x31 - x3*x62 + x30*x9 + x60*x7);
T x65 = x17*x41;
T x66 = x41*(x11*x62 - x13*x31 + x16*x31 - x4*x60 - x6*x61 + x6*x63);
T x67 = x41*(-x10*x11*x33 + x10*x19*x4 + x14*x31 - x19*x6*x7 + x3*x33*x6 - x31*x5);
f[0] = -t - x - y - z + 1;
f[1] = x;
f[2] = y;
f[3] = z;
f[4] = t;
measure_value = x36;
gx[0] = -x42 - x44 - x45;
gy[0] = -x50 - x51 - x54 - x55;
gz[0] = -x56 - x57 - x58 - x59;
gt[0] = -x64 - x65 - x66 - x67;
gx[1] = x45;
gy[1] = x51;
gz[1] = x59;
gt[1] = x65;
gx[2] = x44;
gy[2] = x55;
gz[2] = x57;
gt[2] = x67;
gx[3] = 0;
gy[3] = x50;
gz[3] = x58;
gt[3] = x64;
gx[4] = x42;
gy[4] = x54;
gz[4] = x56;
gt[4] = x66;
			}

		};
	}
}

#endif // UTOPIA_TPL_FE_Pentatope5_4_IMPL_hpp
