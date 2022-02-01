#ifndef UTOPIA_TPL_HYPERELASTICITY_NeohookeanOgden_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeohookeanOgden_3_IMPL_hpp

#include "utopia_hyperelasticity_NeohookeanOgden.hpp"

namespace utopia {
	namespace kernels {

		/** 
		 * Specialization of NeohookeanOgden for dimension 3 
		 */
		template<typename T>
		class NeohookeanOgden<T, 3> {
		public:

			UTOPIA_FUNCTION static void hessian(const T mu,
				const T lmbda,
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *grad_trial,
				const T dx,
				T *UTOPIA_RESTRICT bf)
			{
				using namespace utopia::device;
				// Automatically generated
				T x0 = f[4]*f[8];
				T x1 = f[5]*f[7];
				T x2 = x0 - x1;
				T x3 = f[5]*f[6];
				T x4 = f[3]*f[7];
				T x5 = f[3]*f[8];
				T x6 = f[4]*f[6];
				T x7 = f[0]*x0 - f[0]*x1 + f[1]*x3 - f[1]*x5 + f[2]*x4 - f[2]*x6;
				T x8 = pow(x7, -2);
				T x9 = lmbda*x8;
				T x10 = mu*x8;
				T x11 = -x0 + x1;
				T x12 = x11*x2;
				T x13 = log(x7);
				T x14 = x13*x9;
				T x15 = x3 - x5;
				T x16 = x2*x9;
				T x17 = x15*x16;
				T x18 = x10*x11;
				T x19 = x15*x9;
				T x20 = x11*x13;
				T x21 = x4 - x6;
				T x22 = x16*x21;
				T x23 = x21*x9;
				T x24 = -x3 + x5;
				T x25 = x10*x24;
				T x26 = x13*x24;
				T x27 = x19*x21;
				T x28 = -x4 + x6;
				T x29 = x10*x28;
				T x30 = x13*x28;
				T x31 = f[2]*f[7];
				T x32 = f[1]*f[8];
				T x33 = x31 - x32;
				T x34 = x16*x33;
				T x35 = -x31 + x32;
				T x36 = x10*x35;
				T x37 = x13*x35;
				T x38 = 1.0/x7;
				T x39 = mu*x38;
				T x40 = f[8]*x39;
				T x41 = lmbda*x13*x38;
				T x42 = f[8]*x41;
				T x43 = x19*x33 + x40 - x42;
				T x44 = f[7]*x39;
				T x45 = f[7]*x41;
				T x46 = x23*x33 - x44 + x45;
				T x47 = f[0]*f[8];
				T x48 = f[2]*f[6];
				T x49 = x47 - x48;
				T x50 = x19*x49;
				T x51 = -x47 + x48;
				T x52 = x10*x51;
				T x53 = x13*x51;
				T x54 = x16*x49 - x40 + x42;
				T x55 = f[6]*x39;
				T x56 = f[6]*x41;
				T x57 = x23*x49 + x55 - x56;
				T x58 = f[1]*f[6];
				T x59 = f[0]*f[7];
				T x60 = x58 - x59;
				T x61 = x23*x60;
				T x62 = -x58 + x59;
				T x63 = x10*x62;
				T x64 = x13*x62;
				T x65 = x16*x60 + x44 - x45;
				T x66 = x19*x60 - x55 + x56;
				T x67 = f[1]*f[5];
				T x68 = f[2]*f[4];
				T x69 = x67 - x68;
				T x70 = x16*x69;
				T x71 = -x67 + x68;
				T x72 = x10*x71;
				T x73 = x13*x71;
				T x74 = f[5]*x39;
				T x75 = f[5]*x41;
				T x76 = x19*x69 - x74 + x75;
				T x77 = f[4]*x39;
				T x78 = f[4]*x41;
				T x79 = x23*x69 + x77 - x78;
				T x80 = f[2]*f[3];
				T x81 = f[0]*f[5];
				T x82 = x80 - x81;
				T x83 = x19*x82;
				T x84 = -x80 + x81;
				T x85 = x10*x84;
				T x86 = x13*x84;
				T x87 = x16*x82 + x74 - x75;
				T x88 = f[3]*x39;
				T x89 = f[3]*x41;
				T x90 = x23*x82 - x88 + x89;
				T x91 = f[0]*f[4];
				T x92 = f[1]*f[3];
				T x93 = x91 - x92;
				T x94 = x23*x93;
				T x95 = -x91 + x92;
				T x96 = x10*x95;
				T x97 = x13*x95;
				T x98 = x16*x93 - x77 + x78;
				T x99 = x19*x93 + x88 - x89;
				T x100 = x33*x9;
				T x101 = x49*x9;
				T x102 = x60*x9;
				T x103 = x100*x49;
				T x104 = x100*x60;
				T x105 = x101*x60;
				T x106 = x100*x69;
				T x107 = f[2]*x39;
				T x108 = f[2]*x41;
				T x109 = x101*x69 + x107 - x108;
				T x110 = f[1]*x39;
				T x111 = f[1]*x41;
				T x112 = x102*x69 - x110 + x111;
				T x113 = x101*x82;
				T x114 = x100*x82 - x107 + x108;
				T x115 = f[0]*x39;
				T x116 = f[0]*x41;
				T x117 = x102*x82 + x115 - x116;
				T x118 = x102*x93;
				T x119 = x100*x93 + x110 - x111;
				T x120 = x101*x93 - x115 + x116;
				T x121 = x69*x9;
				T x122 = x82*x9;
				T x123 = x14*x93;
				T x124 = x121*x82;
				T x125 = x121*x93;
				T x126 = x122*x93;
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(mu - x10*x12 + x12*x14 + pow(x2, 2)*x9) + grad_test[1]*(-x15*x18 + x17 + x19*x20) + grad_test[2]*(-x18*x21 + x20*x23 + x22)) + grad_trial[1]*(grad_test[0]*(x16*x26 + x17 - x2*x25) + grad_test[1]*(mu + pow(x15, 2)*x9 - x15*x25 + x19*x26) + grad_test[2]*(-x21*x25 + x23*x26 + x27)) + grad_trial[2]*(grad_test[0]*(x16*x30 - x2*x29 + x22) + grad_test[1]*(-x15*x29 + x19*x30 + x27) + grad_test[2]*(mu + pow(x21, 2)*x9 - x21*x29 + x23*x30)));
				bf[1] += dx*(grad_trial[0]*(grad_test[0]*(x16*x37 - x2*x36 + x34) + grad_test[1]*(-x15*x36 + x19*x37 + x43) + grad_test[2]*(-x21*x36 + x23*x37 + x46)) + grad_trial[1]*(grad_test[0]*(x16*x53 - x2*x52 + x54) + grad_test[1]*(-x15*x52 + x19*x53 + x50) + grad_test[2]*(-x21*x52 + x23*x53 + x57)) + grad_trial[2]*(grad_test[0]*(x16*x64 - x2*x63 + x65) + grad_test[1]*(-x15*x63 + x19*x64 + x66) + grad_test[2]*(-x21*x63 + x23*x64 + x61)));
				bf[2] += dx*(grad_trial[0]*(grad_test[0]*(x16*x73 - x2*x72 + x70) + grad_test[1]*(-x15*x72 + x19*x73 + x76) + grad_test[2]*(-x21*x72 + x23*x73 + x79)) + grad_trial[1]*(grad_test[0]*(x16*x86 - x2*x85 + x87) + grad_test[1]*(-x15*x85 + x19*x86 + x83) + grad_test[2]*(-x21*x85 + x23*x86 + x90)) + grad_trial[2]*(grad_test[0]*(x16*x97 - x2*x96 + x98) + grad_test[1]*(-x15*x96 + x19*x97 + x99) + grad_test[2]*(-x21*x96 + x23*x97 + x94)));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x100*x20 - x18*x33 + x34) + grad_test[1]*(x101*x20 - x18*x49 + x54) + grad_test[2]*(x102*x20 - x18*x60 + x65)) + grad_trial[1]*(grad_test[0]*(x100*x26 - x25*x33 + x43) + grad_test[1]*(x101*x26 - x25*x49 + x50) + grad_test[2]*(x102*x26 - x25*x60 + x66)) + grad_trial[2]*(grad_test[0]*(x100*x30 - x29*x33 + x46) + grad_test[1]*(x101*x30 - x29*x49 + x57) + grad_test[2]*(x102*x30 - x29*x60 + x61)));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(mu + x100*x37 + pow(x33, 2)*x9 - x33*x36) + grad_test[1]*(x101*x37 + x103 - x36*x49) + grad_test[2]*(x102*x37 + x104 - x36*x60)) + grad_trial[1]*(grad_test[0]*(x100*x53 + x103 - x33*x52) + grad_test[1]*(mu + x101*x53 + pow(x49, 2)*x9 - x49*x52) + grad_test[2]*(x102*x53 + x105 - x52*x60)) + grad_trial[2]*(grad_test[0]*(x100*x64 + x104 - x33*x63) + grad_test[1]*(x101*x64 + x105 - x49*x63) + grad_test[2]*(mu + x102*x64 + pow(x60, 2)*x9 - x60*x63)));
				bf[5] += dx*(grad_trial[0]*(grad_test[0]*(x100*x73 + x106 - x33*x72) + grad_test[1]*(x101*x73 + x109 - x49*x72) + grad_test[2]*(x102*x73 + x112 - x60*x72)) + grad_trial[1]*(grad_test[0]*(x100*x86 + x114 - x33*x85) + grad_test[1]*(x101*x86 + x113 - x49*x85) + grad_test[2]*(x102*x86 + x117 - x60*x85)) + grad_trial[2]*(grad_test[0]*(x100*x97 + x119 - x33*x96) + grad_test[1]*(x101*x97 + x120 - x49*x96) + grad_test[2]*(x102*x97 + x118 - x60*x96)));
				bf[6] += dx*(grad_trial[0]*(grad_test[0]*(x121*x20 - x18*x69 + x70) + grad_test[1]*(x122*x20 - x18*x82 + x87) + grad_test[2]*(x11*x123 - x18*x93 + x98)) + grad_trial[1]*(grad_test[0]*(x121*x26 - x25*x69 + x76) + grad_test[1]*(x122*x26 - x25*x82 + x83) + grad_test[2]*(x123*x24 - x25*x93 + x99)) + grad_trial[2]*(grad_test[0]*(x121*x30 - x29*x69 + x79) + grad_test[1]*(x122*x30 - x29*x82 + x90) + grad_test[2]*(x123*x28 - x29*x93 + x94)));
				bf[7] += dx*(grad_trial[0]*(grad_test[0]*(x106 + x121*x37 - x36*x69) + grad_test[1]*(x114 + x122*x37 - x36*x82) + grad_test[2]*(x119 + x123*x35 - x36*x93)) + grad_trial[1]*(grad_test[0]*(x109 + x121*x53 - x52*x69) + grad_test[1]*(x113 + x122*x53 - x52*x82) + grad_test[2]*(x120 + x123*x51 - x52*x93)) + grad_trial[2]*(grad_test[0]*(x112 + x121*x64 - x63*x69) + grad_test[1]*(x117 + x122*x64 - x63*x82) + grad_test[2]*(x118 + x123*x62 - x63*x93)));
				bf[8] += dx*(grad_trial[0]*(grad_test[0]*(mu + x121*x73 + pow(x69, 2)*x9 - x69*x72) + grad_test[1]*(x122*x73 + x124 - x72*x82) + grad_test[2]*(x123*x71 + x125 - x72*x93)) + grad_trial[1]*(grad_test[0]*(x121*x86 + x124 - x69*x85) + grad_test[1]*(mu + x122*x86 + pow(x82, 2)*x9 - x82*x85) + grad_test[2]*(x123*x84 + x126 - x85*x93)) + grad_trial[2]*(grad_test[0]*(x121*x97 + x125 - x69*x96) + grad_test[1]*(x122*x97 + x126 - x82*x96) + grad_test[2]*(mu + x123*x95 + x9*pow(x93, 2) - x93*x96)));
			}

			UTOPIA_FUNCTION static void gradient(const T mu,
				const T lmbda,
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf)
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = f[4]*f[8];
				T x1 = f[5]*f[7];
				T x2 = x0 - x1;
				T x3 = f[5]*f[6];
				T x4 = f[3]*f[7];
				T x5 = f[3]*f[8];
				T x6 = f[4]*f[6];
				T x7 = f[0]*x0 - f[0]*x1 + f[1]*x3 - f[1]*x5 + f[2]*x4 - f[2]*x6;
				T x8 = 1.0/x7;
				T x9 = mu*x8;
				T x10 = lmbda*x8*log(x7);
				T x11 = x3 - x5;
				T x12 = x4 - x6;
				T x13 = -f[1]*f[8] + f[2]*f[7];
				T x14 = f[0]*f[8] - f[2]*f[6];
				T x15 = -f[0]*f[7] + f[1]*f[6];
				T x16 = f[1]*f[5] - f[2]*f[4];
				T x17 = -f[0]*f[5] + f[2]*f[3];
				T x18 = f[0]*f[4] - f[1]*f[3];
				lf[0] += dx*(grad_test[0]*(f[0]*mu + x10*x2 - x2*x9) + grad_test[1]*(f[1]*mu + x10*x11 - x11*x9) + grad_test[2]*(f[2]*mu + x10*x12 - x12*x9));
				lf[1] += dx*(grad_test[0]*(f[3]*mu + x10*x13 - x13*x9) + grad_test[1]*(f[4]*mu + x10*x14 - x14*x9) + grad_test[2]*(f[5]*mu + x10*x15 - x15*x9));
				lf[2] += dx*(grad_test[0]*(f[6]*mu + x10*x16 - x16*x9) + grad_test[1]*(f[7]*mu + x10*x17 - x17*x9) + grad_test[2]*(f[8]*mu + x10*x18 - x18*x9));
			}

			UTOPIA_FUNCTION static void value(const T mu,
				const T lmbda,
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				)
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = log(f[0]*f[4]*f[8] - f[0]*f[5]*f[7] - f[1]*f[3]*f[8] + f[1]*f[5]*f[6] + f[2]*f[3]*f[7] - f[2]*f[4]*f[6]);
				e += dx*((1.0/2.0)*lmbda*pow(x0, 2) - mu*x0 + (1.0/2.0)*mu*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) - 3));
			}

			UTOPIA_FUNCTION void eval(const T mu,
				const T lmbda,
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *grad_trial,
				const T dx,
				T &e,
				T *UTOPIA_RESTRICT lf,
				T *UTOPIA_RESTRICT bf)
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = f[4]*f[8];
				T x1 = f[5]*f[6];
				T x2 = f[3]*f[7];
				T x3 = f[5]*f[7];
				T x4 = f[3]*f[8];
				T x5 = f[4]*f[6];
				T x6 = f[0]*x0 - f[0]*x3 + f[1]*x1 - f[1]*x4 + f[2]*x2 - f[2]*x5;
				T x7 = log(x6);
				T x8 = f[0]*mu;
				T x9 = x0 - x3;
				T x10 = 1.0/x6;
				T x11 = mu*x10;
				T x12 = x7*x9;
				T x13 = lmbda*x10;
				T x14 = f[1]*mu;
				T x15 = x1 - x4;
				T x16 = x13*x7;
				T x17 = f[2]*mu;
				T x18 = x2 - x5;
				T x19 = f[3]*mu;
				T x20 = f[2]*f[7];
				T x21 = f[1]*f[8];
				T x22 = x20 - x21;
				T x23 = f[4]*mu;
				T x24 = f[0]*f[8];
				T x25 = f[2]*f[6];
				T x26 = x24 - x25;
				T x27 = f[5]*mu;
				T x28 = f[1]*f[6];
				T x29 = f[0]*f[7];
				T x30 = x28 - x29;
				T x31 = f[6]*mu;
				T x32 = f[1]*f[5];
				T x33 = f[2]*f[4];
				T x34 = x32 - x33;
				T x35 = f[7]*mu;
				T x36 = f[2]*f[3];
				T x37 = f[0]*f[5];
				T x38 = x36 - x37;
				T x39 = f[8]*mu;
				T x40 = f[0]*f[4];
				T x41 = f[1]*f[3];
				T x42 = x40 - x41;
				T x43 = pow(x6, -2);
				T x44 = lmbda*x43;
				T x45 = -x0 + x3;
				T x46 = mu*x43;
				T x47 = x45*x46;
				T x48 = x12*x44;
				T x49 = x44*x9;
				T x50 = x15*x49;
				T x51 = x15*x44;
				T x52 = x45*x7;
				T x53 = x18*x49;
				T x54 = x18*x44;
				T x55 = -x1 + x4;
				T x56 = x46*x55;
				T x57 = x55*x7;
				T x58 = x18*x51;
				T x59 = -x2 + x5;
				T x60 = x46*x59;
				T x61 = x59*x7;
				T x62 = x22*x49;
				T x63 = -x20 + x21;
				T x64 = x46*x63;
				T x65 = x63*x7;
				T x66 = x10*x39;
				T x67 = f[8]*x16;
				T x68 = x22*x51 + x66 - x67;
				T x69 = x10*x35;
				T x70 = f[7]*x16;
				T x71 = x22*x54 - x69 + x70;
				T x72 = x26*x51;
				T x73 = -x24 + x25;
				T x74 = x46*x73;
				T x75 = x7*x73;
				T x76 = x26*x49 - x66 + x67;
				T x77 = x10*x31;
				T x78 = f[6]*x16;
				T x79 = x26*x54 + x77 - x78;
				T x80 = x30*x54;
				T x81 = -x28 + x29;
				T x82 = x46*x81;
				T x83 = x7*x81;
				T x84 = x30*x49 + x69 - x70;
				T x85 = x30*x51 - x77 + x78;
				T x86 = x34*x49;
				T x87 = -x32 + x33;
				T x88 = x46*x87;
				T x89 = x7*x87;
				T x90 = x10*x27;
				T x91 = f[5]*x16;
				T x92 = x34*x51 - x90 + x91;
				T x93 = x10*x23;
				T x94 = f[4]*x16;
				T x95 = x34*x54 + x93 - x94;
				T x96 = x38*x51;
				T x97 = -x36 + x37;
				T x98 = x46*x97;
				T x99 = x7*x97;
				T x100 = x38*x49 + x90 - x91;
				T x101 = x10*x19;
				T x102 = f[3]*x16;
				T x103 = -x101 + x102 + x38*x54;
				T x104 = x42*x54;
				T x105 = -x40 + x41;
				T x106 = x105*x46;
				T x107 = x105*x7;
				T x108 = x42*x49 - x93 + x94;
				T x109 = x101 - x102 + x42*x51;
				T x110 = x22*x44;
				T x111 = x26*x44;
				T x112 = x30*x44;
				T x113 = x110*x26;
				T x114 = x110*x30;
				T x115 = x111*x30;
				T x116 = x110*x34;
				T x117 = x10*x17;
				T x118 = f[2]*x16;
				T x119 = x111*x34 + x117 - x118;
				T x120 = x10*x14;
				T x121 = f[1]*x16;
				T x122 = x112*x34 - x120 + x121;
				T x123 = x111*x38;
				T x124 = x110*x38 - x117 + x118;
				T x125 = x10*x8;
				T x126 = f[0]*x16;
				T x127 = x112*x38 + x125 - x126;
				T x128 = x112*x42;
				T x129 = x110*x42 + x120 - x121;
				T x130 = x111*x42 - x125 + x126;
				T x131 = x34*x44;
				T x132 = x38*x44;
				T x133 = x42*x44;
				T x134 = x131*x38;
				T x135 = x131*x42;
				T x136 = x132*x42;
				e += dx*((1.0/2.0)*lmbda*pow(x7, 2) - mu*x7 + (1.0/2.0)*mu*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) - 3));
				lf[0] += dx*(grad_test[0]*(-x11*x9 + x12*x13 + x8) + grad_test[1]*(-x11*x15 + x14 + x15*x16) + grad_test[2]*(-x11*x18 + x16*x18 + x17));
				lf[1] += dx*(grad_test[0]*(-x11*x22 + x16*x22 + x19) + grad_test[1]*(-x11*x26 + x16*x26 + x23) + grad_test[2]*(-x11*x30 + x16*x30 + x27));
				lf[2] += dx*(grad_test[0]*(-x11*x34 + x16*x34 + x31) + grad_test[1]*(-x11*x38 + x16*x38 + x35) + grad_test[2]*(-x11*x42 + x16*x42 + x39));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(mu + x44*pow(x9, 2) + x45*x48 - x47*x9) + grad_test[1]*(-x15*x47 + x50 + x51*x52) + grad_test[2]*(-x18*x47 + x52*x54 + x53)) + grad_trial[1]*(grad_test[0]*(x48*x55 + x50 - x56*x9) + grad_test[1]*(mu + pow(x15, 2)*x44 - x15*x56 + x51*x57) + grad_test[2]*(-x18*x56 + x54*x57 + x58)) + grad_trial[2]*(grad_test[0]*(x48*x59 + x53 - x60*x9) + grad_test[1]*(-x15*x60 + x51*x61 + x58) + grad_test[2]*(mu + pow(x18, 2)*x44 - x18*x60 + x54*x61)));
				bf[1] += dx*(grad_trial[0]*(grad_test[0]*(x48*x63 + x62 - x64*x9) + grad_test[1]*(-x15*x64 + x51*x65 + x68) + grad_test[2]*(-x18*x64 + x54*x65 + x71)) + grad_trial[1]*(grad_test[0]*(x48*x73 - x74*x9 + x76) + grad_test[1]*(-x15*x74 + x51*x75 + x72) + grad_test[2]*(-x18*x74 + x54*x75 + x79)) + grad_trial[2]*(grad_test[0]*(x48*x81 - x82*x9 + x84) + grad_test[1]*(-x15*x82 + x51*x83 + x85) + grad_test[2]*(-x18*x82 + x54*x83 + x80)));
				bf[2] += dx*(grad_trial[0]*(grad_test[0]*(x48*x87 + x86 - x88*x9) + grad_test[1]*(-x15*x88 + x51*x89 + x92) + grad_test[2]*(-x18*x88 + x54*x89 + x95)) + grad_trial[1]*(grad_test[0]*(x100 + x48*x97 - x9*x98) + grad_test[1]*(-x15*x98 + x51*x99 + x96) + grad_test[2]*(x103 - x18*x98 + x54*x99)) + grad_trial[2]*(grad_test[0]*(x105*x48 - x106*x9 + x108) + grad_test[1]*(-x106*x15 + x107*x51 + x109) + grad_test[2]*(x104 - x106*x18 + x107*x54)));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x110*x52 - x22*x47 + x62) + grad_test[1]*(x111*x52 - x26*x47 + x76) + grad_test[2]*(x112*x52 - x30*x47 + x84)) + grad_trial[1]*(grad_test[0]*(x110*x57 - x22*x56 + x68) + grad_test[1]*(x111*x57 - x26*x56 + x72) + grad_test[2]*(x112*x57 - x30*x56 + x85)) + grad_trial[2]*(grad_test[0]*(x110*x61 - x22*x60 + x71) + grad_test[1]*(x111*x61 - x26*x60 + x79) + grad_test[2]*(x112*x61 - x30*x60 + x80)));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(mu + x110*x65 + pow(x22, 2)*x44 - x22*x64) + grad_test[1]*(x111*x65 + x113 - x26*x64) + grad_test[2]*(x112*x65 + x114 - x30*x64)) + grad_trial[1]*(grad_test[0]*(x110*x75 + x113 - x22*x74) + grad_test[1]*(mu + x111*x75 + pow(x26, 2)*x44 - x26*x74) + grad_test[2]*(x112*x75 + x115 - x30*x74)) + grad_trial[2]*(grad_test[0]*(x110*x83 + x114 - x22*x82) + grad_test[1]*(x111*x83 + x115 - x26*x82) + grad_test[2]*(mu + x112*x83 + pow(x30, 2)*x44 - x30*x82)));
				bf[5] += dx*(grad_trial[0]*(grad_test[0]*(x110*x89 + x116 - x22*x88) + grad_test[1]*(x111*x89 + x119 - x26*x88) + grad_test[2]*(x112*x89 + x122 - x30*x88)) + grad_trial[1]*(grad_test[0]*(x110*x99 + x124 - x22*x98) + grad_test[1]*(x111*x99 + x123 - x26*x98) + grad_test[2]*(x112*x99 + x127 - x30*x98)) + grad_trial[2]*(grad_test[0]*(-x106*x22 + x107*x110 + x129) + grad_test[1]*(-x106*x26 + x107*x111 + x130) + grad_test[2]*(-x106*x30 + x107*x112 + x128)));
				bf[6] += dx*(grad_trial[0]*(grad_test[0]*(x131*x52 - x34*x47 + x86) + grad_test[1]*(x100 + x132*x52 - x38*x47) + grad_test[2]*(x108 + x133*x52 - x42*x47)) + grad_trial[1]*(grad_test[0]*(x131*x57 - x34*x56 + x92) + grad_test[1]*(x132*x57 - x38*x56 + x96) + grad_test[2]*(x109 + x133*x57 - x42*x56)) + grad_trial[2]*(grad_test[0]*(x131*x61 - x34*x60 + x95) + grad_test[1]*(x103 + x132*x61 - x38*x60) + grad_test[2]*(x104 + x133*x61 - x42*x60)));
				bf[7] += dx*(grad_trial[0]*(grad_test[0]*(x116 + x131*x65 - x34*x64) + grad_test[1]*(x124 + x132*x65 - x38*x64) + grad_test[2]*(x129 + x133*x65 - x42*x64)) + grad_trial[1]*(grad_test[0]*(x119 + x131*x75 - x34*x74) + grad_test[1]*(x123 + x132*x75 - x38*x74) + grad_test[2]*(x130 + x133*x75 - x42*x74)) + grad_trial[2]*(grad_test[0]*(x122 + x131*x83 - x34*x82) + grad_test[1]*(x127 + x132*x83 - x38*x82) + grad_test[2]*(x128 + x133*x83 - x42*x82)));
				bf[8] += dx*(grad_trial[0]*(grad_test[0]*(mu + x131*x89 + pow(x34, 2)*x44 - x34*x88) + grad_test[1]*(x132*x89 + x134 - x38*x88) + grad_test[2]*(x133*x89 + x135 - x42*x88)) + grad_trial[1]*(grad_test[0]*(x131*x99 + x134 - x34*x98) + grad_test[1]*(mu + x132*x99 + pow(x38, 2)*x44 - x38*x98) + grad_test[2]*(x133*x99 + x136 - x42*x98)) + grad_trial[2]*(grad_test[0]*(-x106*x34 + x107*x131 + x135) + grad_test[1]*(-x106*x38 + x107*x132 + x136) + grad_test[2]*(mu - x106*x42 + x107*x133 + pow(x42, 2)*x44)));
			}

		};;
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_NeohookeanOgden_3_IMPL_hpp
