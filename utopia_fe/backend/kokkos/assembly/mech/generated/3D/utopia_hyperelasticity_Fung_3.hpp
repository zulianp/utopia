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
			class Params : public Configurable {
			public:
				void read(Input &in) override
				{
					in.get("a", a);
					in.get("b", b);
					in.get("c", c);
				}

				T a;
				T b;
				T c;
			};

			Fung(const Params &params)
			{
				a = params.a;
				b = params.b;
				c = params.c;
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
				T x0 = pow(f[0], 2);
				T x1 = pow(f[1], 2);
				T x2 = pow(f[2], 2);
				T x3 = pow(f[3], 2);
				T x4 = pow(f[4], 2);
				T x5 = pow(f[5], 2);
				T x6 = pow(f[6], 2);
				T x7 = pow(f[7], 2);
				T x8 = pow(f[8], 2);
				T x9 = b*exp(c*(x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 - 3) - 1);
				T x10 = 2*pow(c, 2)*x9;
				T x11 = f[0]*x10;
				T x12 = f[1]*x11;
				T x13 = f[2]*grad_test[2];
				T x14 = a + c*x9;
				T x15 = f[1]*x10;
				T x16 = grad_test[0]*x11;
				T x17 = grad_test[1]*x15;
				T x18 = f[3]*x16;
				T x19 = x10*x13;
				T x20 = f[4]*x17;
				T x21 = f[5]*x19;
				T x22 = f[6]*x16;
				T x23 = f[7]*x17;
				T x24 = f[8]*x19;
				T x25 = f[4]*grad_test[1];
				T x26 = f[5]*grad_test[2];
				T x27 = f[3]*grad_test[0];
				T x28 = f[2]*x10;
				T x29 = f[3]*x10;
				T x30 = f[4]*x10;
				T x31 = f[5]*x10;
				T x32 = f[6]*x10;
				T x33 = x27*x32;
				T x34 = f[7]*x10;
				T x35 = x25*x34;
				T x36 = f[8]*x10;
				T x37 = x26*x36;
				T x38 = f[7]*grad_test[1];
				T x39 = f[8]*grad_test[2];
				T x40 = f[6]*grad_test[0];
				T x41 = grad_test[0]*x32;
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x0*x10 + x14) + grad_test[1]*x12 + x11*x13) + grad_trial[1]*(grad_test[0]*x12 + grad_test[1]*(x1*x10 + x14) + x13*x15) + grad_trial[2]*(f[2]*x16 + f[2]*x17 + grad_test[2]*(x10*x2 + x14)));
				bf[1] += dx*(grad_trial[0]*(f[3]*x17 + f[3]*x19 + x18) + grad_trial[1]*(f[4]*x16 + f[4]*x19 + x20) + grad_trial[2]*(f[5]*x16 + f[5]*x17 + x21));
				bf[2] += dx*(grad_trial[0]*(f[6]*x17 + f[6]*x19 + x22) + grad_trial[1]*(f[7]*x16 + f[7]*x19 + x23) + grad_trial[2]*(f[8]*x16 + f[8]*x17 + x24));
				bf[3] += dx*(grad_trial[0]*(x11*x25 + x11*x26 + x18) + grad_trial[1]*(x15*x26 + x15*x27 + x20) + grad_trial[2]*(x21 + x25*x28 + x27*x28));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(x10*x3 + x14) + x25*x29 + x26*x29) + grad_trial[1]*(grad_test[1]*(x10*x4 + x14) + x26*x30 + x27*x30) + grad_trial[2]*(grad_test[2]*(x10*x5 + x14) + x25*x31 + x27*x31));
				bf[5] += dx*(grad_trial[0]*(x25*x32 + x26*x32 + x33) + grad_trial[1]*(x26*x34 + x27*x34 + x35) + grad_trial[2]*(x25*x36 + x27*x36 + x37));
				bf[6] += dx*(grad_trial[0]*(x11*x38 + x11*x39 + x22) + grad_trial[1]*(x15*x39 + x15*x40 + x23) + grad_trial[2]*(x24 + x28*x38 + x28*x40));
				bf[7] += dx*(grad_trial[0]*(x29*x38 + x29*x39 + x33) + grad_trial[1]*(x30*x39 + x30*x40 + x35) + grad_trial[2]*(x31*x38 + x31*x40 + x37));
				bf[8] += dx*(grad_trial[0]*(grad_test[0]*(x10*x6 + x14) + x32*x38 + x32*x39) + grad_trial[1]*(f[7]*x41 + grad_test[1]*(x10*x7 + x14) + x34*x39) + grad_trial[2]*(f[8]*grad_test[1]*x34 + f[8]*x41 + grad_test[2]*(x10*x8 + x14)));
			}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T *UTOPIA_RESTRICT grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = b*c*exp(c*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) - 3) - 1);
				lf[0] += dx*(grad_test[0]*(a*f[0] + f[0]*x0) + grad_test[1]*(a*f[1] + f[1]*x0) + grad_test[2]*(a*f[2] + f[2]*x0));
				lf[1] += dx*(grad_test[0]*(a*f[3] + f[3]*x0) + grad_test[1]*(a*f[4] + f[4]*x0) + grad_test[2]*(a*f[5] + f[5]*x0));
				lf[2] += dx*(grad_test[0]*(a*f[6] + f[6]*x0) + grad_test[1]*(a*f[7] + f[7]*x0) + grad_test[2]*(a*f[8] + f[8]*x0));
			}

			UTOPIA_FUNCTION void value(
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) - 3;
				e += dx*((1.0/2.0)*a*x0 + (1.0/2.0)*b*exp(c*x0 - 1));
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
				T x0 = pow(f[0], 2);
				T x1 = pow(f[1], 2);
				T x2 = pow(f[2], 2);
				T x3 = pow(f[3], 2);
				T x4 = pow(f[4], 2);
				T x5 = pow(f[5], 2);
				T x6 = pow(f[6], 2);
				T x7 = pow(f[7], 2);
				T x8 = pow(f[8], 2);
				T x9 = x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 - 3;
				T x10 = b*exp(c*x9 - 1);
				T x11 = c*x10;
				T x12 = 2*pow(c, 2)*x10;
				T x13 = f[0]*x12;
				T x14 = f[1]*x13;
				T x15 = f[2]*grad_test[2];
				T x16 = a + x11;
				T x17 = f[1]*x12;
				T x18 = grad_test[0]*x13;
				T x19 = grad_test[1]*x17;
				T x20 = f[3]*x18;
				T x21 = x12*x15;
				T x22 = f[4]*x19;
				T x23 = f[5]*x21;
				T x24 = f[6]*x18;
				T x25 = f[7]*x19;
				T x26 = f[8]*x21;
				T x27 = f[4]*grad_test[1];
				T x28 = f[5]*grad_test[2];
				T x29 = f[3]*grad_test[0];
				T x30 = f[2]*x12;
				T x31 = f[3]*x12;
				T x32 = f[4]*x12;
				T x33 = f[5]*x12;
				T x34 = f[6]*x12;
				T x35 = x29*x34;
				T x36 = f[7]*x12;
				T x37 = x27*x36;
				T x38 = f[8]*x12;
				T x39 = x28*x38;
				T x40 = f[7]*grad_test[1];
				T x41 = f[8]*grad_test[2];
				T x42 = f[6]*grad_test[0];
				T x43 = grad_test[0]*x34;
				e += dx*((1.0/2.0)*a*x9 + (1.0/2.0)*x10);
				lf[0] += dx*(grad_test[0]*(a*f[0] + f[0]*x11) + grad_test[1]*(a*f[1] + f[1]*x11) + grad_test[2]*(a*f[2] + f[2]*x11));
				lf[1] += dx*(grad_test[0]*(a*f[3] + f[3]*x11) + grad_test[1]*(a*f[4] + f[4]*x11) + grad_test[2]*(a*f[5] + f[5]*x11));
				lf[2] += dx*(grad_test[0]*(a*f[6] + f[6]*x11) + grad_test[1]*(a*f[7] + f[7]*x11) + grad_test[2]*(a*f[8] + f[8]*x11));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x0*x12 + x16) + grad_test[1]*x14 + x13*x15) + grad_trial[1]*(grad_test[0]*x14 + grad_test[1]*(x1*x12 + x16) + x15*x17) + grad_trial[2]*(f[2]*x18 + f[2]*x19 + grad_test[2]*(x12*x2 + x16)));
				bf[1] += dx*(grad_trial[0]*(f[3]*x19 + f[3]*x21 + x20) + grad_trial[1]*(f[4]*x18 + f[4]*x21 + x22) + grad_trial[2]*(f[5]*x18 + f[5]*x19 + x23));
				bf[2] += dx*(grad_trial[0]*(f[6]*x19 + f[6]*x21 + x24) + grad_trial[1]*(f[7]*x18 + f[7]*x21 + x25) + grad_trial[2]*(f[8]*x18 + f[8]*x19 + x26));
				bf[3] += dx*(grad_trial[0]*(x13*x27 + x13*x28 + x20) + grad_trial[1]*(x17*x28 + x17*x29 + x22) + grad_trial[2]*(x23 + x27*x30 + x29*x30));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(x12*x3 + x16) + x27*x31 + x28*x31) + grad_trial[1]*(grad_test[1]*(x12*x4 + x16) + x28*x32 + x29*x32) + grad_trial[2]*(grad_test[2]*(x12*x5 + x16) + x27*x33 + x29*x33));
				bf[5] += dx*(grad_trial[0]*(x27*x34 + x28*x34 + x35) + grad_trial[1]*(x28*x36 + x29*x36 + x37) + grad_trial[2]*(x27*x38 + x29*x38 + x39));
				bf[6] += dx*(grad_trial[0]*(x13*x40 + x13*x41 + x24) + grad_trial[1]*(x17*x41 + x17*x42 + x25) + grad_trial[2]*(x26 + x30*x40 + x30*x42));
				bf[7] += dx*(grad_trial[0]*(x31*x40 + x31*x41 + x35) + grad_trial[1]*(x32*x41 + x32*x42 + x37) + grad_trial[2]*(x33*x40 + x33*x42 + x39));
				bf[8] += dx*(grad_trial[0]*(grad_test[0]*(x12*x6 + x16) + x34*x40 + x34*x41) + grad_trial[1]*(f[7]*x43 + grad_test[1]*(x12*x7 + x16) + x36*x41) + grad_trial[2]*(f[8]*grad_test[1]*x36 + f[8]*x43 + grad_test[2]*(x12*x8 + x16)));
			}

			T a;
			T b;
			T c;

		};
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_Fung_3_IMPL_hpp
