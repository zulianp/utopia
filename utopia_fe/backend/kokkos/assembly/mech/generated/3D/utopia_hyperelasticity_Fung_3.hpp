#ifndef UTOPIA_TPL_HYPERELASTICITY_Fung_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_Fung_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_Fung.hpp"

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of Fung for dimension 3
		 */
		template<typename T>
		class Fung<T, 3> {
		public:
			static constexpr int Dim = 3;

			class Params : public Configurable {
			public:
				void read(Input &in) override
				{
					in.get("a", a);
					in.get("b", b);
					in.get("c", c);
					in.get("k", k);
				}

				T a{1.0};
				T b{1.0};
				T c{1.0};
				T k{1};
			};

			Fung(const Params &params)
			{
				a = params.a;
				b = params.b;
				c = params.c;
				k = params.k;
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
				T x0 = f[3]*f[8];
				T x1 = f[5]*f[6];
				T x2 = -2*x0 + 2*x1;
				T x3 = f[4]*f[8];
				T x4 = f[5]*f[7];
				T x5 = (1.0/2.0)*k;
				T x6 = x5*(x3 - x4);
				T x7 = pow(b, 2);
				T x8 = (1.0/2.0)*pow(f[0], 2);
				T x9 = (1.0/2.0)*pow(f[1], 2);
				T x10 = (1.0/2.0)*pow(f[2], 2);
				T x11 = (1.0/2.0)*pow(f[3], 2);
				T x12 = (1.0/2.0)*pow(f[4], 2);
				T x13 = (1.0/2.0)*pow(f[5], 2);
				T x14 = (1.0/2.0)*pow(f[6], 2);
				T x15 = (1.0/2.0)*pow(f[7], 2);
				T x16 = (1.0/2.0)*pow(f[8], 2);
				T x17 = c*exp(b*(x10 + x11 + x12 + x13 + x14 + x15 + x16 + x8 + x9 - 3.0/2.0));
				T x18 = (1.0/2.0)*x17;
				T x19 = x18*x7;
				T x20 = f[0]*x19;
				T x21 = f[1]*x20;
				T x22 = f[3]*f[7];
				T x23 = f[4]*f[6];
				T x24 = 2*x22 - 2*x23;
				T x25 = f[2]*x20;
				T x26 = 2*x3 - 2*x4;
				T x27 = x17*x7;
				T x28 = (1.0/2.0)*a + b*x18;
				T x29 = x5*(-x0 + x1);
				T x30 = f[1]*x19;
				T x31 = f[2]*x30;
				T x32 = x5*(x22 - x23);
				T x33 = f[2]*f[7];
				T x34 = f[1]*f[8];
				T x35 = x5*(x33 - x34);
				T x36 = f[3]*x20;
				T x37 = k*(f[0]*x3 - f[0]*x4 - f[1]*x0 + f[1]*x1 + f[2]*x22 - f[2]*x23 - 1);
				T x38 = f[7]*x37;
				T x39 = f[2]*f[3];
				T x40 = x19*x39 + x38;
				T x41 = f[8]*x37;
				T x42 = f[1]*f[3];
				T x43 = x19*x42 - x41;
				T x44 = f[0]*f[8];
				T x45 = f[2]*f[6];
				T x46 = x5*(x44 - x45);
				T x47 = f[4]*x30;
				T x48 = f[0]*f[4];
				T x49 = x19*x48 + x41;
				T x50 = f[6]*x37;
				T x51 = f[2]*f[4];
				T x52 = x19*x51 - x50;
				T x53 = f[1]*f[6];
				T x54 = f[0]*f[7];
				T x55 = x5*(x53 - x54);
				T x56 = f[2]*x19;
				T x57 = f[5]*x56;
				T x58 = f[1]*f[5];
				T x59 = x19*x58 + x50;
				T x60 = f[0]*f[5];
				T x61 = x19*x60 - x38;
				T x62 = x5*(-x51 + x58);
				T x63 = f[6]*x20;
				T x64 = f[5]*x37;
				T x65 = x19*x53 + x64;
				T x66 = f[4]*x37;
				T x67 = x19*x45 - x66;
				T x68 = x5*(x39 - x60);
				T x69 = f[7]*x30;
				T x70 = f[3]*x37;
				T x71 = x19*x33 + x70;
				T x72 = x19*x54 - x64;
				T x73 = x5*(-x42 + x48);
				T x74 = f[8]*x56;
				T x75 = x19*x44 + x66;
				T x76 = x19*x34 - x70;
				T x77 = 2*x33 - 2*x34;
				T x78 = 2*x44 - 2*x45;
				T x79 = 2*x53 - 2*x54;
				T x80 = f[3]*x19;
				T x81 = f[4]*x80;
				T x82 = f[5]*x80;
				T x83 = f[4]*x19;
				T x84 = f[5]*x83;
				T x85 = f[6]*x80;
				T x86 = f[1]*x37;
				T x87 = x1*x19 + x86;
				T x88 = f[2]*x37;
				T x89 = x19*x23 - x88;
				T x90 = f[7]*x83;
				T x91 = x19*x22 + x88;
				T x92 = f[0]*x37;
				T x93 = x19*x4 - x92;
				T x94 = f[8]*x19;
				T x95 = f[5]*x94;
				T x96 = x19*x3 + x92;
				T x97 = x0*x19 - x86;
				T x98 = -2*x51 + 2*x58;
				T x99 = -2*x42 + 2*x48;
				T x100 = 2*x39 - 2*x60;
				T x101 = f[6]*f[7]*x19;
				T x102 = f[6]*x94;
				T x103 = f[7]*x94;
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x26*x6 + x27*x8 + x28) + grad_test[1]*(x2*x6 + x21) + grad_test[2]*(x24*x6 + x25)) + grad_trial[1]*(grad_test[0]*(x21 + x26*x29) + grad_test[1]*(x2*x29 + x27*x9 + x28) + grad_test[2]*(x24*x29 + x31)) + grad_trial[2]*(grad_test[0]*(x25 + x26*x32) + grad_test[1]*(x2*x32 + x31) + grad_test[2]*(x10*x27 + x24*x32 + x28)));
				bf[1] += dx*(grad_trial[0]*(grad_test[0]*(x26*x35 + x36) + grad_test[1]*(x2*x35 + x43) + grad_test[2]*(x24*x35 + x40)) + grad_trial[1]*(grad_test[0]*(x26*x46 + x49) + grad_test[1]*(x2*x46 + x47) + grad_test[2]*(x24*x46 + x52)) + grad_trial[2]*(grad_test[0]*(x26*x55 + x61) + grad_test[1]*(x2*x55 + x59) + grad_test[2]*(x24*x55 + x57)));
				bf[2] += dx*(grad_trial[0]*(grad_test[0]*(x26*x62 + x63) + grad_test[1]*(x2*x62 + x65) + grad_test[2]*(x24*x62 + x67)) + grad_trial[1]*(grad_test[0]*(x26*x68 + x72) + grad_test[1]*(x2*x68 + x69) + grad_test[2]*(x24*x68 + x71)) + grad_trial[2]*(grad_test[0]*(x26*x73 + x75) + grad_test[1]*(x2*x73 + x76) + grad_test[2]*(x24*x73 + x74)));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x36 + x6*x77) + grad_test[1]*(x49 + x6*x78) + grad_test[2]*(x6*x79 + x61)) + grad_trial[1]*(grad_test[0]*(x29*x77 + x43) + grad_test[1]*(x29*x78 + x47) + grad_test[2]*(x29*x79 + x59)) + grad_trial[2]*(grad_test[0]*(x32*x77 + x40) + grad_test[1]*(x32*x78 + x52) + grad_test[2]*(x32*x79 + x57)));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(x11*x27 + x28 + x35*x77) + grad_test[1]*(x35*x78 + x81) + grad_test[2]*(x35*x79 + x82)) + grad_trial[1]*(grad_test[0]*(x46*x77 + x81) + grad_test[1]*(x12*x27 + x28 + x46*x78) + grad_test[2]*(x46*x79 + x84)) + grad_trial[2]*(grad_test[0]*(x55*x77 + x82) + grad_test[1]*(x55*x78 + x84) + grad_test[2]*(x13*x27 + x28 + x55*x79)));
				bf[5] += dx*(grad_trial[0]*(grad_test[0]*(x62*x77 + x85) + grad_test[1]*(x62*x78 + x89) + grad_test[2]*(x62*x79 + x87)) + grad_trial[1]*(grad_test[0]*(x68*x77 + x91) + grad_test[1]*(x68*x78 + x90) + grad_test[2]*(x68*x79 + x93)) + grad_trial[2]*(grad_test[0]*(x73*x77 + x97) + grad_test[1]*(x73*x78 + x96) + grad_test[2]*(x73*x79 + x95)));
				bf[6] += dx*(grad_trial[0]*(grad_test[0]*(x6*x98 + x63) + grad_test[1]*(x100*x6 + x72) + grad_test[2]*(x6*x99 + x75)) + grad_trial[1]*(grad_test[0]*(x29*x98 + x65) + grad_test[1]*(x100*x29 + x69) + grad_test[2]*(x29*x99 + x76)) + grad_trial[2]*(grad_test[0]*(x32*x98 + x67) + grad_test[1]*(x100*x32 + x71) + grad_test[2]*(x32*x99 + x74)));
				bf[7] += dx*(grad_trial[0]*(grad_test[0]*(x35*x98 + x85) + grad_test[1]*(x100*x35 + x91) + grad_test[2]*(x35*x99 + x97)) + grad_trial[1]*(grad_test[0]*(x46*x98 + x89) + grad_test[1]*(x100*x46 + x90) + grad_test[2]*(x46*x99 + x96)) + grad_trial[2]*(grad_test[0]*(x55*x98 + x87) + grad_test[1]*(x100*x55 + x93) + grad_test[2]*(x55*x99 + x95)));
				bf[8] += dx*(grad_trial[0]*(grad_test[0]*(x14*x27 + x28 + x62*x98) + grad_test[1]*(x100*x62 + x101) + grad_test[2]*(x102 + x62*x99)) + grad_trial[1]*(grad_test[0]*(x101 + x68*x98) + grad_test[1]*(x100*x68 + x15*x27 + x28) + grad_test[2]*(x103 + x68*x99)) + grad_trial[2]*(grad_test[0]*(x102 + x73*x98) + grad_test[1]*(x100*x73 + x103) + grad_test[2]*(x16*x27 + x28 + x73*x99)));
			}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T *UTOPIA_RESTRICT grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = (1.0/2.0)*a;
				T x1 = f[4]*f[8];
				T x2 = f[5]*f[7];
				T x3 = f[5]*f[6];
				T x4 = f[3]*f[7];
				T x5 = f[3]*f[8];
				T x6 = f[4]*f[6];
				T x7 = (1.0/2.0)*k*(f[0]*x1 - f[0]*x2 + f[1]*x3 - f[1]*x5 + f[2]*x4 - f[2]*x6 - 1);
				T x8 = (1.0/2.0)*b*c*exp(b*((1.0/2.0)*pow(f[0], 2) + (1.0/2.0)*pow(f[1], 2) + (1.0/2.0)*pow(f[2], 2) + (1.0/2.0)*pow(f[3], 2) + (1.0/2.0)*pow(f[4], 2) + (1.0/2.0)*pow(f[5], 2) + (1.0/2.0)*pow(f[6], 2) + (1.0/2.0)*pow(f[7], 2) + (1.0/2.0)*pow(f[8], 2) - 3.0/2.0));
				T x9 = 2*f[8];
				T x10 = 2*f[2];
				T x11 = 2*f[0];
				T x12 = 2*f[1];
				lf[0] += dx*(grad_test[0]*(f[0]*x0 + f[0]*x8 + x7*(2*x1 - 2*x2)) + grad_test[1]*(f[1]*x0 + f[1]*x8 + x7*(2*x3 - 2*x5)) + grad_test[2]*(f[2]*x0 + f[2]*x8 + x7*(2*x4 - 2*x6)));
				lf[1] += dx*(grad_test[0]*(f[3]*x0 + f[3]*x8 + x7*(-f[1]*x9 + f[7]*x10)) + grad_test[1]*(f[4]*x0 + f[4]*x8 + x7*(f[0]*x9 - f[6]*x10)) + grad_test[2]*(f[5]*x0 + f[5]*x8 + x7*(f[6]*x12 - f[7]*x11)));
				lf[2] += dx*(grad_test[0]*(f[6]*x0 + f[6]*x8 + x7*(-f[4]*x10 + f[5]*x12)) + grad_test[1]*(f[7]*x0 + f[7]*x8 + x7*(f[3]*x10 - f[5]*x11)) + grad_test[2]*(f[8]*x0 + f[8]*x8 + x7*(-f[3]*x12 + f[4]*x11)));
			}

			UTOPIA_FUNCTION void value(
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = (1.0/2.0)*pow(f[0], 2) + (1.0/2.0)*pow(f[1], 2) + (1.0/2.0)*pow(f[2], 2) + (1.0/2.0)*pow(f[3], 2) + (1.0/2.0)*pow(f[4], 2) + (1.0/2.0)*pow(f[5], 2) + (1.0/2.0)*pow(f[6], 2) + (1.0/2.0)*pow(f[7], 2) + (1.0/2.0)*pow(f[8], 2) - 3.0/2.0;
				e += dx*((1.0/2.0)*a*x0 + (1.0/2.0)*c*(exp(b*x0) - 1) + (1.0/2.0)*k*pow(f[0]*f[4]*f[8] - f[0]*f[5]*f[7] - f[1]*f[3]*f[8] + f[1]*f[5]*f[6] + f[2]*f[3]*f[7] - f[2]*f[4]*f[6] - 1, 2));
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
				T x6 = f[0]*x0 - f[0]*x3 + f[1]*x1 - f[1]*x4 + f[2]*x2 - f[2]*x5 - 1;
				T x7 = (1.0/2.0)*k;
				T x8 = (1.0/2.0)*pow(f[0], 2);
				T x9 = (1.0/2.0)*pow(f[1], 2);
				T x10 = (1.0/2.0)*pow(f[2], 2);
				T x11 = (1.0/2.0)*pow(f[3], 2);
				T x12 = (1.0/2.0)*pow(f[4], 2);
				T x13 = (1.0/2.0)*pow(f[5], 2);
				T x14 = (1.0/2.0)*pow(f[6], 2);
				T x15 = (1.0/2.0)*pow(f[7], 2);
				T x16 = (1.0/2.0)*pow(f[8], 2);
				T x17 = x10 + x11 + x12 + x13 + x14 + x15 + x16 + x8 + x9 - 3.0/2.0;
				T x18 = (1.0/2.0)*a;
				T x19 = exp(b*x17);
				T x20 = (1.0/2.0)*c;
				T x21 = 2*x0 - 2*x3;
				T x22 = x6*x7;
				T x23 = x19*x20;
				T x24 = b*x23;
				T x25 = 2*x1 - 2*x4;
				T x26 = 2*x2 - 2*x5;
				T x27 = f[1]*f[8];
				T x28 = f[2]*f[7];
				T x29 = -2*x27 + 2*x28;
				T x30 = f[0]*f[8];
				T x31 = f[2]*f[6];
				T x32 = 2*x30 - 2*x31;
				T x33 = f[0]*f[7];
				T x34 = f[1]*f[6];
				T x35 = -2*x33 + 2*x34;
				T x36 = f[1]*f[5];
				T x37 = f[2]*f[4];
				T x38 = 2*x36 - 2*x37;
				T x39 = f[0]*f[5];
				T x40 = f[2]*f[3];
				T x41 = -2*x39 + 2*x40;
				T x42 = f[0]*f[4];
				T x43 = f[1]*f[3];
				T x44 = 2*x42 - 2*x43;
				T x45 = x7*(x0 - x3);
				T x46 = pow(b, 2);
				T x47 = x23*x46;
				T x48 = f[0]*x47;
				T x49 = f[1]*x48;
				T x50 = f[2]*x48;
				T x51 = c*x19*x46;
				T x52 = x18 + x24;
				T x53 = x7*(x1 - x4);
				T x54 = f[1]*x47;
				T x55 = f[2]*x54;
				T x56 = x7*(x2 - x5);
				T x57 = x7*(-x27 + x28);
				T x58 = f[3]*x48;
				T x59 = k*x6;
				T x60 = f[7]*x59;
				T x61 = x40*x47 + x60;
				T x62 = f[8]*x59;
				T x63 = x43*x47 - x62;
				T x64 = x7*(x30 - x31);
				T x65 = f[4]*x54;
				T x66 = x42*x47 + x62;
				T x67 = f[6]*x59;
				T x68 = x37*x47 - x67;
				T x69 = x7*(-x33 + x34);
				T x70 = f[2]*x47;
				T x71 = f[5]*x70;
				T x72 = x36*x47 + x67;
				T x73 = x39*x47 - x60;
				T x74 = x7*(x36 - x37);
				T x75 = f[6]*x48;
				T x76 = f[5]*x59;
				T x77 = x34*x47 + x76;
				T x78 = f[4]*x59;
				T x79 = x31*x47 - x78;
				T x80 = x7*(-x39 + x40);
				T x81 = f[7]*x54;
				T x82 = f[3]*x59;
				T x83 = x28*x47 + x82;
				T x84 = x33*x47 - x76;
				T x85 = x7*(x42 - x43);
				T x86 = f[8]*x70;
				T x87 = x30*x47 + x78;
				T x88 = x27*x47 - x82;
				T x89 = f[3]*x47;
				T x90 = f[4]*x89;
				T x91 = f[5]*x89;
				T x92 = f[4]*x47;
				T x93 = f[5]*x92;
				T x94 = f[6]*x89;
				T x95 = f[1]*x59;
				T x96 = x1*x47 + x95;
				T x97 = f[2]*x59;
				T x98 = x47*x5 - x97;
				T x99 = f[7]*x92;
				T x100 = x2*x47 + x97;
				T x101 = f[0]*x59;
				T x102 = -x101 + x3*x47;
				T x103 = f[8]*x47;
				T x104 = f[5]*x103;
				T x105 = x0*x47 + x101;
				T x106 = x4*x47 - x95;
				T x107 = f[6]*f[7]*x47;
				T x108 = f[6]*x103;
				T x109 = f[7]*x103;
				e += dx*(x17*x18 + x20*(x19 - 1) + pow(x6, 2)*x7);
				lf[0] += dx*(grad_test[0]*(f[0]*x18 + f[0]*x24 + x21*x22) + grad_test[1]*(f[1]*x18 + f[1]*x24 + x22*x25) + grad_test[2]*(f[2]*x18 + f[2]*x24 + x22*x26));
				lf[1] += dx*(grad_test[0]*(f[3]*x18 + f[3]*x24 + x22*x29) + grad_test[1]*(f[4]*x18 + f[4]*x24 + x22*x32) + grad_test[2]*(f[5]*x18 + f[5]*x24 + x22*x35));
				lf[2] += dx*(grad_test[0]*(f[6]*x18 + f[6]*x24 + x22*x38) + grad_test[1]*(f[7]*x18 + f[7]*x24 + x22*x41) + grad_test[2]*(f[8]*x18 + f[8]*x24 + x22*x44));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x21*x45 + x51*x8 + x52) + grad_test[1]*(x25*x45 + x49) + grad_test[2]*(x26*x45 + x50)) + grad_trial[1]*(grad_test[0]*(x21*x53 + x49) + grad_test[1]*(x25*x53 + x51*x9 + x52) + grad_test[2]*(x26*x53 + x55)) + grad_trial[2]*(grad_test[0]*(x21*x56 + x50) + grad_test[1]*(x25*x56 + x55) + grad_test[2]*(x10*x51 + x26*x56 + x52)));
				bf[1] += dx*(grad_trial[0]*(grad_test[0]*(x21*x57 + x58) + grad_test[1]*(x25*x57 + x63) + grad_test[2]*(x26*x57 + x61)) + grad_trial[1]*(grad_test[0]*(x21*x64 + x66) + grad_test[1]*(x25*x64 + x65) + grad_test[2]*(x26*x64 + x68)) + grad_trial[2]*(grad_test[0]*(x21*x69 + x73) + grad_test[1]*(x25*x69 + x72) + grad_test[2]*(x26*x69 + x71)));
				bf[2] += dx*(grad_trial[0]*(grad_test[0]*(x21*x74 + x75) + grad_test[1]*(x25*x74 + x77) + grad_test[2]*(x26*x74 + x79)) + grad_trial[1]*(grad_test[0]*(x21*x80 + x84) + grad_test[1]*(x25*x80 + x81) + grad_test[2]*(x26*x80 + x83)) + grad_trial[2]*(grad_test[0]*(x21*x85 + x87) + grad_test[1]*(x25*x85 + x88) + grad_test[2]*(x26*x85 + x86)));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x29*x45 + x58) + grad_test[1]*(x32*x45 + x66) + grad_test[2]*(x35*x45 + x73)) + grad_trial[1]*(grad_test[0]*(x29*x53 + x63) + grad_test[1]*(x32*x53 + x65) + grad_test[2]*(x35*x53 + x72)) + grad_trial[2]*(grad_test[0]*(x29*x56 + x61) + grad_test[1]*(x32*x56 + x68) + grad_test[2]*(x35*x56 + x71)));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(x11*x51 + x29*x57 + x52) + grad_test[1]*(x32*x57 + x90) + grad_test[2]*(x35*x57 + x91)) + grad_trial[1]*(grad_test[0]*(x29*x64 + x90) + grad_test[1]*(x12*x51 + x32*x64 + x52) + grad_test[2]*(x35*x64 + x93)) + grad_trial[2]*(grad_test[0]*(x29*x69 + x91) + grad_test[1]*(x32*x69 + x93) + grad_test[2]*(x13*x51 + x35*x69 + x52)));
				bf[5] += dx*(grad_trial[0]*(grad_test[0]*(x29*x74 + x94) + grad_test[1]*(x32*x74 + x98) + grad_test[2]*(x35*x74 + x96)) + grad_trial[1]*(grad_test[0]*(x100 + x29*x80) + grad_test[1]*(x32*x80 + x99) + grad_test[2]*(x102 + x35*x80)) + grad_trial[2]*(grad_test[0]*(x106 + x29*x85) + grad_test[1]*(x105 + x32*x85) + grad_test[2]*(x104 + x35*x85)));
				bf[6] += dx*(grad_trial[0]*(grad_test[0]*(x38*x45 + x75) + grad_test[1]*(x41*x45 + x84) + grad_test[2]*(x44*x45 + x87)) + grad_trial[1]*(grad_test[0]*(x38*x53 + x77) + grad_test[1]*(x41*x53 + x81) + grad_test[2]*(x44*x53 + x88)) + grad_trial[2]*(grad_test[0]*(x38*x56 + x79) + grad_test[1]*(x41*x56 + x83) + grad_test[2]*(x44*x56 + x86)));
				bf[7] += dx*(grad_trial[0]*(grad_test[0]*(x38*x57 + x94) + grad_test[1]*(x100 + x41*x57) + grad_test[2]*(x106 + x44*x57)) + grad_trial[1]*(grad_test[0]*(x38*x64 + x98) + grad_test[1]*(x41*x64 + x99) + grad_test[2]*(x105 + x44*x64)) + grad_trial[2]*(grad_test[0]*(x38*x69 + x96) + grad_test[1]*(x102 + x41*x69) + grad_test[2]*(x104 + x44*x69)));
				bf[8] += dx*(grad_trial[0]*(grad_test[0]*(x14*x51 + x38*x74 + x52) + grad_test[1]*(x107 + x41*x74) + grad_test[2]*(x108 + x44*x74)) + grad_trial[1]*(grad_test[0]*(x107 + x38*x80) + grad_test[1]*(x15*x51 + x41*x80 + x52) + grad_test[2]*(x109 + x44*x80)) + grad_trial[2]*(grad_test[0]*(x108 + x38*x85) + grad_test[1]*(x109 + x41*x85) + grad_test[2]*(x16*x51 + x44*x85 + x52)));
			}

			T a{1.0};
			T b{1.0};
			T c{1.0};
			T k{1};

		};
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_Fung_3_IMPL_hpp
