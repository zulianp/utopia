#ifndef UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeohookeanBower.hpp"

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of NeohookeanBower for dimension 3
		 */
		template<typename T>
		class NeohookeanBower<T, 3> {
		public:
			class Params : public Configurable {
			public:
				void read(Input &in) override
				{
					in.get("mu", mu);
					in.get("lambda", lambda);
				}

				T mu;
				T lambda;
			};

			NeohookeanBower(const Params &params)
			{
				mu = params.mu;
				lambda = params.lambda;
			}

			UTOPIA_FUNCTION void hessian(
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *grad_trial,
				const T dx,
				T *UTOPIA_RESTRICT bf) const
			{
				using namespace utopia::device;
				// Automatically generated
				T x0 = f[4]*f[8];
				T x1 = f[5]*f[7];
				T x2 = 2*x0 - 2*x1;
				T x3 = (1.0/2.0)*lambda;
				T x4 = x3*(x0 - x1);
				T x5 = f[5]*f[6];
				T x6 = f[3]*f[7];
				T x7 = f[3]*f[8];
				T x8 = f[4]*f[6];
				T x9 = f[0]*x0 - f[0]*x1 + f[1]*x5 - f[1]*x7 + f[2]*x6 - f[2]*x8;
				T x10 = 2/pow(x9, 2.0/3.0);
				T x11 = -2.0/3.0*x0 + (2.0/3.0)*x1;
				T x12 = pow(x9, -5.0/3.0);
				T x13 = f[0]*x12;
				T x14 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2);
				T x15 = x14/pow(x9, 8.0/3.0);
				T x16 = x15*(-5.0/3.0*x0 + (5.0/3.0)*x1);
				T x17 = (1.0/2.0)*mu;
				T x18 = 2*x5 - 2*x7;
				T x19 = -2.0/3.0*x5 + (2.0/3.0)*x7;
				T x20 = 2*x13;
				T x21 = f[1]*x12;
				T x22 = 2*x11;
				T x23 = x19*x20 + x21*x22;
				T x24 = 2*x6 - 2*x8;
				T x25 = -2.0/3.0*x6 + (2.0/3.0)*x8;
				T x26 = f[2]*x12;
				T x27 = x20*x25 + x22*x26;
				T x28 = x3*(x5 - x7);
				T x29 = x15*(-5.0/3.0*x5 + (5.0/3.0)*x7);
				T x30 = 2*x21;
				T x31 = 2*x19;
				T x32 = x25*x30 + x26*x31;
				T x33 = x3*(x6 - x8);
				T x34 = x15*(-5.0/3.0*x6 + (5.0/3.0)*x8);
				T x35 = f[2]*f[7];
				T x36 = f[1]*f[8];
				T x37 = x3*(x35 - x36);
				T x38 = x15*(-5.0/3.0*x35 + (5.0/3.0)*x36);
				T x39 = -2.0/3.0*x35 + (2.0/3.0)*x36;
				T x40 = f[3]*x12;
				T x41 = x20*x39 + x22*x40;
				T x42 = lambda*(x9 - 1);
				T x43 = f[7]*x42;
				T x44 = 2*x26;
				T x45 = 2*x25;
				T x46 = (2.0/3.0)*x14;
				T x47 = x12*x46;
				T x48 = f[7]*x47;
				T x49 = x39*x44 + x40*x45 - x48;
				T x50 = f[8]*x42;
				T x51 = -x50;
				T x52 = f[8]*x12;
				T x53 = x46*x52;
				T x54 = x30*x39 + x31*x40 + x53;
				T x55 = f[0]*f[8];
				T x56 = f[2]*f[6];
				T x57 = x3*(x55 - x56);
				T x58 = x15*(-5.0/3.0*x55 + (5.0/3.0)*x56);
				T x59 = -2.0/3.0*x55 + (2.0/3.0)*x56;
				T x60 = f[4]*x12;
				T x61 = x30*x59 + x31*x60;
				T x62 = x20*x59 + x22*x60 - x53;
				T x63 = f[6]*x42;
				T x64 = -x63;
				T x65 = f[6]*x47;
				T x66 = x44*x59 + x45*x60 + x65;
				T x67 = f[1]*f[6];
				T x68 = f[0]*f[7];
				T x69 = x3*(x67 - x68);
				T x70 = x15*(-5.0/3.0*x67 + (5.0/3.0)*x68);
				T x71 = -2.0/3.0*x67 + (2.0/3.0)*x68;
				T x72 = f[5]*x12;
				T x73 = x44*x71 + x45*x72;
				T x74 = x30*x71 + x31*x72 - x65;
				T x75 = -x43;
				T x76 = x20*x71 + x22*x72 + x48;
				T x77 = f[1]*f[5];
				T x78 = f[2]*f[4];
				T x79 = x3*(x77 - x78);
				T x80 = x15*(-5.0/3.0*x77 + (5.0/3.0)*x78);
				T x81 = -2.0/3.0*x77 + (2.0/3.0)*x78;
				T x82 = f[6]*x12;
				T x83 = x20*x81 + x22*x82;
				T x84 = f[5]*x42;
				T x85 = f[5]*x47;
				T x86 = x30*x81 + x31*x82 - x85;
				T x87 = f[4]*x42;
				T x88 = -x87;
				T x89 = f[4]*x47;
				T x90 = x44*x81 + x45*x82 + x89;
				T x91 = f[2]*f[3];
				T x92 = f[0]*f[5];
				T x93 = x3*(x91 - x92);
				T x94 = x15*(-5.0/3.0*x91 + (5.0/3.0)*x92);
				T x95 = -2.0/3.0*x91 + (2.0/3.0)*x92;
				T x96 = f[7]*x12;
				T x97 = x30*x95 + x31*x96;
				T x98 = f[3]*x42;
				T x99 = x40*x46;
				T x100 = x44*x95 + x45*x96 - x99;
				T x101 = -x84;
				T x102 = x20*x95 + x22*x96 + x85;
				T x103 = f[0]*f[4];
				T x104 = f[1]*f[3];
				T x105 = x3*(x103 - x104);
				T x106 = x15*(-5.0/3.0*x103 + (5.0/3.0)*x104);
				T x107 = -2.0/3.0*x103 + (2.0/3.0)*x104;
				T x108 = x107*x44 + x45*x52;
				T x109 = x107*x20 + x22*x52 - x89;
				T x110 = -x98;
				T x111 = x107*x30 + x31*x52 + x99;
				T x112 = 2*x35 - 2*x36;
				T x113 = 2*x55 - 2*x56;
				T x114 = 2*x67 - 2*x68;
				T x115 = 2*x40;
				T x116 = 2*x39;
				T x117 = x115*x59 + x116*x60;
				T x118 = x115*x71 + x116*x72;
				T x119 = 2*x60;
				T x120 = 2*x59;
				T x121 = x119*x71 + x120*x72;
				T x122 = x115*x81 + x116*x82;
				T x123 = f[1]*x42;
				T x124 = 2*x72;
				T x125 = 2*x71;
				T x126 = x21*x46;
				T x127 = x124*x81 + x125*x82 - x126;
				T x128 = f[2]*x42;
				T x129 = -x128;
				T x130 = x26*x46;
				T x131 = x119*x81 + x120*x82 + x130;
				T x132 = x119*x95 + x120*x96;
				T x133 = x115*x95 + x116*x96 - x130;
				T x134 = f[0]*x42;
				T x135 = -x134;
				T x136 = x13*x46;
				T x137 = x124*x95 + x125*x96 + x136;
				T x138 = x107*x124 + x125*x52;
				T x139 = x107*x119 + x120*x52 - x136;
				T x140 = -x123;
				T x141 = x107*x115 + x116*x52 + x126;
				T x142 = 2*x77 - 2*x78;
				T x143 = 2*x103 - 2*x104;
				T x144 = 2*x91 - 2*x92;
				T x145 = 2*x82;
				T x146 = 2*x81;
				T x147 = x145*x95 + x146*x96;
				T x148 = x107*x145 + x146*x52;
				T x149 = 2*x107*x96 + 2*x52*x95;
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x17*(x10 + 4*x11*x13 + x11*x16) + x2*x4) + grad_test[1]*(x17*(x16*x19 + x23) + x18*x4) + grad_test[2]*(x17*(x16*x25 + x27) + x24*x4)) + grad_trial[1]*(grad_test[0]*(x17*(x11*x29 + x23) + x2*x28) + grad_test[1]*(x17*(x10 + 4*x19*x21 + x19*x29) + x18*x28) + grad_test[2]*(x17*(x25*x29 + x32) + x24*x28)) + grad_trial[2]*(grad_test[0]*(x17*(x11*x34 + x27) + x2*x33) + grad_test[1]*(x17*(x19*x34 + x32) + x18*x33) + grad_test[2]*(x17*(x10 + 4*x25*x26 + x25*x34) + x24*x33)));
				bf[1] += dx*(grad_trial[0]*(grad_test[0]*(x17*(x11*x38 + x41) + x2*x37) + grad_test[1]*(x17*(x19*x38 + x54) + x18*x37 + x51) + grad_test[2]*(x17*(x25*x38 + x49) + x24*x37 + x43)) + grad_trial[1]*(grad_test[0]*(x17*(x11*x58 + x62) + x2*x57 + x50) + grad_test[1]*(x17*(x19*x58 + x61) + x18*x57) + grad_test[2]*(x17*(x25*x58 + x66) + x24*x57 + x64)) + grad_trial[2]*(grad_test[0]*(x17*(x11*x70 + x76) + x2*x69 + x75) + grad_test[1]*(x17*(x19*x70 + x74) + x18*x69 + x63) + grad_test[2]*(x17*(x25*x70 + x73) + x24*x69)));
				bf[2] += dx*(grad_trial[0]*(grad_test[0]*(x17*(x11*x80 + x83) + x2*x79) + grad_test[1]*(x17*(x19*x80 + x86) + x18*x79 + x84) + grad_test[2]*(x17*(x25*x80 + x90) + x24*x79 + x88)) + grad_trial[1]*(grad_test[0]*(x101 + x17*(x102 + x11*x94) + x2*x93) + grad_test[1]*(x17*(x19*x94 + x97) + x18*x93) + grad_test[2]*(x17*(x100 + x25*x94) + x24*x93 + x98)) + grad_trial[2]*(grad_test[0]*(x105*x2 + x17*(x106*x11 + x109) + x87) + grad_test[1]*(x105*x18 + x110 + x17*(x106*x19 + x111)) + grad_test[2]*(x105*x24 + x17*(x106*x25 + x108))));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x112*x4 + x17*(x16*x39 + x41)) + grad_test[1]*(x113*x4 + x17*(x16*x59 + x62) + x50) + grad_test[2]*(x114*x4 + x17*(x16*x71 + x76) + x75)) + grad_trial[1]*(grad_test[0]*(x112*x28 + x17*(x29*x39 + x54) + x51) + grad_test[1]*(x113*x28 + x17*(x29*x59 + x61)) + grad_test[2]*(x114*x28 + x17*(x29*x71 + x74) + x63)) + grad_trial[2]*(grad_test[0]*(x112*x33 + x17*(x34*x39 + x49) + x43) + grad_test[1]*(x113*x33 + x17*(x34*x59 + x66) + x64) + grad_test[2]*(x114*x33 + x17*(x34*x71 + x73))));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(x112*x37 + x17*(x10 + x38*x39 + 4*x39*x40)) + grad_test[1]*(x113*x37 + x17*(x117 + x38*x59)) + grad_test[2]*(x114*x37 + x17*(x118 + x38*x71))) + grad_trial[1]*(grad_test[0]*(x112*x57 + x17*(x117 + x39*x58)) + grad_test[1]*(x113*x57 + x17*(x10 + x58*x59 + 4*x59*x60)) + grad_test[2]*(x114*x57 + x17*(x121 + x58*x71))) + grad_trial[2]*(grad_test[0]*(x112*x69 + x17*(x118 + x39*x70)) + grad_test[1]*(x113*x69 + x17*(x121 + x59*x70)) + grad_test[2]*(x114*x69 + x17*(x10 + x70*x71 + 4*x71*x72))));
				bf[5] += dx*(grad_trial[0]*(grad_test[0]*(x112*x79 + x17*(x122 + x39*x80)) + grad_test[1]*(x113*x79 + x129 + x17*(x131 + x59*x80)) + grad_test[2]*(x114*x79 + x123 + x17*(x127 + x71*x80))) + grad_trial[1]*(grad_test[0]*(x112*x93 + x128 + x17*(x133 + x39*x94)) + grad_test[1]*(x113*x93 + x17*(x132 + x59*x94)) + grad_test[2]*(x114*x93 + x135 + x17*(x137 + x71*x94))) + grad_trial[2]*(grad_test[0]*(x105*x112 + x140 + x17*(x106*x39 + x141)) + grad_test[1]*(x105*x113 + x134 + x17*(x106*x59 + x139)) + grad_test[2]*(x105*x114 + x17*(x106*x71 + x138))));
				bf[6] += dx*(grad_trial[0]*(grad_test[0]*(x142*x4 + x17*(x16*x81 + x83)) + grad_test[1]*(x101 + x144*x4 + x17*(x102 + x16*x95)) + grad_test[2]*(x143*x4 + x17*(x107*x16 + x109) + x87)) + grad_trial[1]*(grad_test[0]*(x142*x28 + x17*(x29*x81 + x86) + x84) + grad_test[1]*(x144*x28 + x17*(x29*x95 + x97)) + grad_test[2]*(x110 + x143*x28 + x17*(x107*x29 + x111))) + grad_trial[2]*(grad_test[0]*(x142*x33 + x17*(x34*x81 + x90) + x88) + grad_test[1]*(x144*x33 + x17*(x100 + x34*x95) + x98) + grad_test[2]*(x143*x33 + x17*(x107*x34 + x108))));
				bf[7] += dx*(grad_trial[0]*(grad_test[0]*(x142*x37 + x17*(x122 + x38*x81)) + grad_test[1]*(x128 + x144*x37 + x17*(x133 + x38*x95)) + grad_test[2]*(x140 + x143*x37 + x17*(x107*x38 + x141))) + grad_trial[1]*(grad_test[0]*(x129 + x142*x57 + x17*(x131 + x58*x81)) + grad_test[1]*(x144*x57 + x17*(x132 + x58*x95)) + grad_test[2]*(x134 + x143*x57 + x17*(x107*x58 + x139))) + grad_trial[2]*(grad_test[0]*(x123 + x142*x69 + x17*(x127 + x70*x81)) + grad_test[1]*(x135 + x144*x69 + x17*(x137 + x70*x95)) + grad_test[2]*(x143*x69 + x17*(x107*x70 + x138))));
				bf[8] += dx*(grad_trial[0]*(grad_test[0]*(x142*x79 + x17*(x10 + x80*x81 + 4*x81*x82)) + grad_test[1]*(x144*x79 + x17*(x147 + x80*x95)) + grad_test[2]*(x143*x79 + x17*(x107*x80 + x148))) + grad_trial[1]*(grad_test[0]*(x142*x93 + x17*(x147 + x81*x94)) + grad_test[1]*(x144*x93 + x17*(x10 + x94*x95 + 4*x95*x96)) + grad_test[2]*(x143*x93 + x17*(x107*x94 + x149))) + grad_trial[2]*(grad_test[0]*(x105*x142 + x17*(x106*x81 + x148)) + grad_test[1]*(x105*x144 + x17*(x106*x95 + x149)) + grad_test[2]*(x105*x143 + x17*(x10 + x106*x107 + 4*x107*x52))));
			}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T *UTOPIA_RESTRICT grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = f[4]*f[8];
				T x1 = f[5]*f[7];
				T x2 = f[5]*f[6];
				T x3 = f[3]*f[7];
				T x4 = f[3]*f[8];
				T x5 = f[4]*f[6];
				T x6 = f[0]*x0 - f[0]*x1 + f[1]*x2 - f[1]*x4 + f[2]*x3 - f[2]*x5;
				T x7 = (1.0/2.0)*lambda*(x6 - 1);
				T x8 = 2/pow(x6, 2.0/3.0);
				T x9 = (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2))/pow(x6, 5.0/3.0);
				T x10 = (1.0/2.0)*mu;
				T x11 = f[1]*f[8];
				T x12 = f[2]*f[7];
				T x13 = f[0]*f[8];
				T x14 = f[2]*f[6];
				T x15 = f[0]*f[7];
				T x16 = f[1]*f[6];
				T x17 = f[1]*f[5];
				T x18 = f[2]*f[4];
				T x19 = f[0]*f[5];
				T x20 = f[2]*f[3];
				T x21 = f[0]*f[4];
				T x22 = f[1]*f[3];
				lf[0] += dx*(grad_test[0]*(x10*(f[0]*x8 + x9*(-2.0/3.0*x0 + (2.0/3.0)*x1)) + x7*(2*x0 - 2*x1)) + grad_test[1]*(x10*(f[1]*x8 + x9*(-2.0/3.0*x2 + (2.0/3.0)*x4)) + x7*(2*x2 - 2*x4)) + grad_test[2]*(x10*(f[2]*x8 + x9*(-2.0/3.0*x3 + (2.0/3.0)*x5)) + x7*(2*x3 - 2*x5)));
				lf[1] += dx*(grad_test[0]*(x10*(f[3]*x8 + x9*((2.0/3.0)*x11 - 2.0/3.0*x12)) + x7*(-2*x11 + 2*x12)) + grad_test[1]*(x10*(f[4]*x8 + x9*(-2.0/3.0*x13 + (2.0/3.0)*x14)) + x7*(2*x13 - 2*x14)) + grad_test[2]*(x10*(f[5]*x8 + x9*((2.0/3.0)*x15 - 2.0/3.0*x16)) + x7*(-2*x15 + 2*x16)));
				lf[2] += dx*(grad_test[0]*(x10*(f[6]*x8 + x9*(-2.0/3.0*x17 + (2.0/3.0)*x18)) + x7*(2*x17 - 2*x18)) + grad_test[1]*(x10*(f[7]*x8 + x9*((2.0/3.0)*x19 - 2.0/3.0*x20)) + x7*(-2*x19 + 2*x20)) + grad_test[2]*(x10*(f[8]*x8 + x9*(-2.0/3.0*x21 + (2.0/3.0)*x22)) + x7*(2*x21 - 2*x22)));
			}

			UTOPIA_FUNCTION void value(
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = f[0]*f[4]*f[8] - f[0]*f[5]*f[7] - f[1]*f[3]*f[8] + f[1]*f[5]*f[6] + f[2]*f[3]*f[7] - f[2]*f[4]*f[6];
				e += dx*((1.0/2.0)*lambda*pow(x0 - 1, 2) + (1.0/2.0)*mu*(-3 + (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2))/pow(x0, 2.0/3.0)));
			}

			UTOPIA_FUNCTION void eval(
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *grad_trial,
				const T dx,
				T &e,
				T *UTOPIA_RESTRICT lf,
				T *UTOPIA_RESTRICT bf) const
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
				T x7 = x6 - 1;
				T x8 = (1.0/2.0)*lambda;
				T x9 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2);
				T x10 = pow(x6, -2.0/3.0);
				T x11 = (1.0/2.0)*mu;
				T x12 = 2*x0 - 2*x3;
				T x13 = x7*x8;
				T x14 = 2*x10;
				T x15 = -2.0/3.0*x0 + (2.0/3.0)*x3;
				T x16 = pow(x6, -5.0/3.0);
				T x17 = x16*x9;
				T x18 = 2*x1 - 2*x4;
				T x19 = -2.0/3.0*x1 + (2.0/3.0)*x4;
				T x20 = 2*x2 - 2*x5;
				T x21 = -2.0/3.0*x2 + (2.0/3.0)*x5;
				T x22 = f[1]*f[8];
				T x23 = f[2]*f[7];
				T x24 = -2*x22 + 2*x23;
				T x25 = (2.0/3.0)*x22 - 2.0/3.0*x23;
				T x26 = f[0]*f[8];
				T x27 = f[2]*f[6];
				T x28 = 2*x26 - 2*x27;
				T x29 = -2.0/3.0*x26 + (2.0/3.0)*x27;
				T x30 = f[0]*f[7];
				T x31 = f[1]*f[6];
				T x32 = -2*x30 + 2*x31;
				T x33 = (2.0/3.0)*x30 - 2.0/3.0*x31;
				T x34 = f[1]*f[5];
				T x35 = f[2]*f[4];
				T x36 = 2*x34 - 2*x35;
				T x37 = -2.0/3.0*x34 + (2.0/3.0)*x35;
				T x38 = f[0]*f[5];
				T x39 = f[2]*f[3];
				T x40 = -2*x38 + 2*x39;
				T x41 = (2.0/3.0)*x38 - 2.0/3.0*x39;
				T x42 = f[0]*f[4];
				T x43 = f[1]*f[3];
				T x44 = 2*x42 - 2*x43;
				T x45 = -2.0/3.0*x42 + (2.0/3.0)*x43;
				T x46 = x8*(x0 - x3);
				T x47 = f[0]*x16;
				T x48 = x9/pow(x6, 8.0/3.0);
				T x49 = x48*(-5.0/3.0*x0 + (5.0/3.0)*x3);
				T x50 = 2*x47;
				T x51 = f[1]*x16;
				T x52 = 2*x15;
				T x53 = x19*x50 + x51*x52;
				T x54 = f[2]*x16;
				T x55 = x21*x50 + x52*x54;
				T x56 = x8*(x1 - x4);
				T x57 = x48*(-5.0/3.0*x1 + (5.0/3.0)*x4);
				T x58 = 2*x51;
				T x59 = 2*x19;
				T x60 = x21*x58 + x54*x59;
				T x61 = x8*(x2 - x5);
				T x62 = x48*(-5.0/3.0*x2 + (5.0/3.0)*x5);
				T x63 = x8*(-x22 + x23);
				T x64 = x48*((5.0/3.0)*x22 - 5.0/3.0*x23);
				T x65 = f[3]*x16;
				T x66 = x25*x50 + x52*x65;
				T x67 = lambda*x7;
				T x68 = f[7]*x67;
				T x69 = 2*x54;
				T x70 = 2*x21;
				T x71 = (2.0/3.0)*x17;
				T x72 = f[7]*x71;
				T x73 = x25*x69 + x65*x70 - x72;
				T x74 = f[8]*x67;
				T x75 = -x74;
				T x76 = f[8]*x71;
				T x77 = x25*x58 + x59*x65 + x76;
				T x78 = x8*(x26 - x27);
				T x79 = x48*(-5.0/3.0*x26 + (5.0/3.0)*x27);
				T x80 = f[4]*x16;
				T x81 = x29*x58 + x59*x80;
				T x82 = x29*x50 + x52*x80 - x76;
				T x83 = f[6]*x67;
				T x84 = -x83;
				T x85 = f[6]*x71;
				T x86 = x29*x69 + x70*x80 + x85;
				T x87 = x8*(-x30 + x31);
				T x88 = x48*((5.0/3.0)*x30 - 5.0/3.0*x31);
				T x89 = f[5]*x16;
				T x90 = x33*x69 + x70*x89;
				T x91 = x33*x58 + x59*x89 - x85;
				T x92 = -x68;
				T x93 = x33*x50 + x52*x89 + x72;
				T x94 = x8*(x34 - x35);
				T x95 = x48*(-5.0/3.0*x34 + (5.0/3.0)*x35);
				T x96 = f[6]*x16;
				T x97 = x37*x50 + x52*x96;
				T x98 = f[5]*x67;
				T x99 = f[5]*x71;
				T x100 = x37*x58 + x59*x96 - x99;
				T x101 = f[4]*x67;
				T x102 = -x101;
				T x103 = f[4]*x71;
				T x104 = x103 + x37*x69 + x70*x96;
				T x105 = x8*(-x38 + x39);
				T x106 = x48*((5.0/3.0)*x38 - 5.0/3.0*x39);
				T x107 = f[7]*x16;
				T x108 = x107*x59 + x41*x58;
				T x109 = f[3]*x67;
				T x110 = f[3]*x71;
				T x111 = x107*x70 - x110 + x41*x69;
				T x112 = -x98;
				T x113 = x107*x52 + x41*x50 + x99;
				T x114 = x8*(x42 - x43);
				T x115 = x48*(-5.0/3.0*x42 + (5.0/3.0)*x43);
				T x116 = f[8]*x16;
				T x117 = x116*x70 + x45*x69;
				T x118 = -x103 + x116*x52 + x45*x50;
				T x119 = -x109;
				T x120 = x110 + x116*x59 + x45*x58;
				T x121 = 2*x65;
				T x122 = 2*x25;
				T x123 = x121*x29 + x122*x80;
				T x124 = x121*x33 + x122*x89;
				T x125 = 2*x80;
				T x126 = 2*x29;
				T x127 = x125*x33 + x126*x89;
				T x128 = x121*x37 + x122*x96;
				T x129 = f[1]*x67;
				T x130 = 2*x89;
				T x131 = 2*x33;
				T x132 = f[1]*x71;
				T x133 = x130*x37 + x131*x96 - x132;
				T x134 = f[2]*x67;
				T x135 = -x134;
				T x136 = f[2]*x71;
				T x137 = x125*x37 + x126*x96 + x136;
				T x138 = x107*x126 + x125*x41;
				T x139 = x107*x122 + x121*x41 - x136;
				T x140 = f[0]*x67;
				T x141 = -x140;
				T x142 = f[0]*x71;
				T x143 = x107*x131 + x130*x41 + x142;
				T x144 = x116*x131 + x130*x45;
				T x145 = x116*x126 + x125*x45 - x142;
				T x146 = -x129;
				T x147 = x116*x122 + x121*x45 + x132;
				T x148 = 2*x96;
				T x149 = 2*x37;
				T x150 = x107*x149 + x148*x41;
				T x151 = x116*x149 + x148*x45;
				T x152 = 2*x107*x45 + 2*x116*x41;
				e += dx*(x11*(x10*x9 - 3) + pow(x7, 2)*x8);
				lf[0] += dx*(grad_test[0]*(x11*(f[0]*x14 + x15*x17) + x12*x13) + grad_test[1]*(x11*(f[1]*x14 + x17*x19) + x13*x18) + grad_test[2]*(x11*(f[2]*x14 + x17*x21) + x13*x20));
				lf[1] += dx*(grad_test[0]*(x11*(f[3]*x14 + x17*x25) + x13*x24) + grad_test[1]*(x11*(f[4]*x14 + x17*x29) + x13*x28) + grad_test[2]*(x11*(f[5]*x14 + x17*x33) + x13*x32));
				lf[2] += dx*(grad_test[0]*(x11*(f[6]*x14 + x17*x37) + x13*x36) + grad_test[1]*(x11*(f[7]*x14 + x17*x41) + x13*x40) + grad_test[2]*(x11*(f[8]*x14 + x17*x45) + x13*x44));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x14 + 4*x15*x47 + x15*x49) + x12*x46) + grad_test[1]*(x11*(x19*x49 + x53) + x18*x46) + grad_test[2]*(x11*(x21*x49 + x55) + x20*x46)) + grad_trial[1]*(grad_test[0]*(x11*(x15*x57 + x53) + x12*x56) + grad_test[1]*(x11*(x14 + 4*x19*x51 + x19*x57) + x18*x56) + grad_test[2]*(x11*(x21*x57 + x60) + x20*x56)) + grad_trial[2]*(grad_test[0]*(x11*(x15*x62 + x55) + x12*x61) + grad_test[1]*(x11*(x19*x62 + x60) + x18*x61) + grad_test[2]*(x11*(x14 + 4*x21*x54 + x21*x62) + x20*x61)));
				bf[1] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x15*x64 + x66) + x12*x63) + grad_test[1]*(x11*(x19*x64 + x77) + x18*x63 + x75) + grad_test[2]*(x11*(x21*x64 + x73) + x20*x63 + x68)) + grad_trial[1]*(grad_test[0]*(x11*(x15*x79 + x82) + x12*x78 + x74) + grad_test[1]*(x11*(x19*x79 + x81) + x18*x78) + grad_test[2]*(x11*(x21*x79 + x86) + x20*x78 + x84)) + grad_trial[2]*(grad_test[0]*(x11*(x15*x88 + x93) + x12*x87 + x92) + grad_test[1]*(x11*(x19*x88 + x91) + x18*x87 + x83) + grad_test[2]*(x11*(x21*x88 + x90) + x20*x87)));
				bf[2] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x15*x95 + x97) + x12*x94) + grad_test[1]*(x11*(x100 + x19*x95) + x18*x94 + x98) + grad_test[2]*(x102 + x11*(x104 + x21*x95) + x20*x94)) + grad_trial[1]*(grad_test[0]*(x105*x12 + x11*(x106*x15 + x113) + x112) + grad_test[1]*(x105*x18 + x11*(x106*x19 + x108)) + grad_test[2]*(x105*x20 + x109 + x11*(x106*x21 + x111))) + grad_trial[2]*(grad_test[0]*(x101 + x11*(x115*x15 + x118) + x114*x12) + grad_test[1]*(x11*(x115*x19 + x120) + x114*x18 + x119) + grad_test[2]*(x11*(x115*x21 + x117) + x114*x20)));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x25*x49 + x66) + x24*x46) + grad_test[1]*(x11*(x29*x49 + x82) + x28*x46 + x74) + grad_test[2]*(x11*(x33*x49 + x93) + x32*x46 + x92)) + grad_trial[1]*(grad_test[0]*(x11*(x25*x57 + x77) + x24*x56 + x75) + grad_test[1]*(x11*(x29*x57 + x81) + x28*x56) + grad_test[2]*(x11*(x33*x57 + x91) + x32*x56 + x83)) + grad_trial[2]*(grad_test[0]*(x11*(x25*x62 + x73) + x24*x61 + x68) + grad_test[1]*(x11*(x29*x62 + x86) + x28*x61 + x84) + grad_test[2]*(x11*(x33*x62 + x90) + x32*x61)));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x14 + x25*x64 + 4*x25*x65) + x24*x63) + grad_test[1]*(x11*(x123 + x29*x64) + x28*x63) + grad_test[2]*(x11*(x124 + x33*x64) + x32*x63)) + grad_trial[1]*(grad_test[0]*(x11*(x123 + x25*x79) + x24*x78) + grad_test[1]*(x11*(x14 + x29*x79 + 4*x29*x80) + x28*x78) + grad_test[2]*(x11*(x127 + x33*x79) + x32*x78)) + grad_trial[2]*(grad_test[0]*(x11*(x124 + x25*x88) + x24*x87) + grad_test[1]*(x11*(x127 + x29*x88) + x28*x87) + grad_test[2]*(x11*(x14 + x33*x88 + 4*x33*x89) + x32*x87)));
				bf[5] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x128 + x25*x95) + x24*x94) + grad_test[1]*(x11*(x137 + x29*x95) + x135 + x28*x94) + grad_test[2]*(x11*(x133 + x33*x95) + x129 + x32*x94)) + grad_trial[1]*(grad_test[0]*(x105*x24 + x11*(x106*x25 + x139) + x134) + grad_test[1]*(x105*x28 + x11*(x106*x29 + x138)) + grad_test[2]*(x105*x32 + x11*(x106*x33 + x143) + x141)) + grad_trial[2]*(grad_test[0]*(x11*(x115*x25 + x147) + x114*x24 + x146) + grad_test[1]*(x11*(x115*x29 + x145) + x114*x28 + x140) + grad_test[2]*(x11*(x115*x33 + x144) + x114*x32)));
				bf[6] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x37*x49 + x97) + x36*x46) + grad_test[1]*(x11*(x113 + x41*x49) + x112 + x40*x46) + grad_test[2]*(x101 + x11*(x118 + x45*x49) + x44*x46)) + grad_trial[1]*(grad_test[0]*(x11*(x100 + x37*x57) + x36*x56 + x98) + grad_test[1]*(x11*(x108 + x41*x57) + x40*x56) + grad_test[2]*(x11*(x120 + x45*x57) + x119 + x44*x56)) + grad_trial[2]*(grad_test[0]*(x102 + x11*(x104 + x37*x62) + x36*x61) + grad_test[1]*(x109 + x11*(x111 + x41*x62) + x40*x61) + grad_test[2]*(x11*(x117 + x45*x62) + x44*x61)));
				bf[7] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x128 + x37*x64) + x36*x63) + grad_test[1]*(x11*(x139 + x41*x64) + x134 + x40*x63) + grad_test[2]*(x11*(x147 + x45*x64) + x146 + x44*x63)) + grad_trial[1]*(grad_test[0]*(x11*(x137 + x37*x79) + x135 + x36*x78) + grad_test[1]*(x11*(x138 + x41*x79) + x40*x78) + grad_test[2]*(x11*(x145 + x45*x79) + x140 + x44*x78)) + grad_trial[2]*(grad_test[0]*(x11*(x133 + x37*x88) + x129 + x36*x87) + grad_test[1]*(x11*(x143 + x41*x88) + x141 + x40*x87) + grad_test[2]*(x11*(x144 + x45*x88) + x44*x87)));
				bf[8] += dx*(grad_trial[0]*(grad_test[0]*(x11*(x14 + x37*x95 + 4*x37*x96) + x36*x94) + grad_test[1]*(x11*(x150 + x41*x95) + x40*x94) + grad_test[2]*(x11*(x151 + x45*x95) + x44*x94)) + grad_trial[1]*(grad_test[0]*(x105*x36 + x11*(x106*x37 + x150)) + grad_test[1]*(x105*x40 + x11*(x106*x41 + 4*x107*x41 + x14)) + grad_test[2]*(x105*x44 + x11*(x106*x45 + x152))) + grad_trial[2]*(grad_test[0]*(x11*(x115*x37 + x151) + x114*x36) + grad_test[1]*(x11*(x115*x41 + x152) + x114*x40) + grad_test[2]*(x11*(x115*x45 + 4*x116*x45 + x14) + x114*x44)));
			}

			T mu;
			T lambda;

		};
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_3_IMPL_hpp
