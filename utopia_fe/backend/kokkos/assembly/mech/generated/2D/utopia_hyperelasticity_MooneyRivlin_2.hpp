#ifndef UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_MooneyRivlin.hpp"

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of MooneyRivlin for dimension 2
		 */
		template<typename T>
		class MooneyRivlin<T, 2> {
		public:
			static constexpr int Dim = 2;

			class Params : public Configurable {
			public:
				void read(Input &in) override
				{
					in.get("C1", C1);
					in.get("C2", C2);
				}

				T C1{1.0};
				T C2{1.0};
			};

			MooneyRivlin(const Params &params)
			{
				C1 = params.C1;
				C2 = params.C2;
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
				T x0 = pow(f[3], 2);
				T x1 = pow(f[0], 2);
				T x2 = pow(f[2], 2);
				T x3 = x1 + x2;
				T x4 = pow(f[1], 2);
				T x5 = x0 + x4;
				T x6 = x3 + x5;
				T x7 = f[0]*f[3];
				T x8 = f[1]*f[2];
				T x9 = x7 - x8;
				T x10 = pow(x9, -5.0/2.0);
				T x11 = (3.0/4.0)*x10*x6;
				T x12 = 2/sqrt(x9);
				T x13 = pow(x9, -3.0/2.0);
				T x14 = 2*x7;
				T x15 = x12 - x13*x14;
				T x16 = x0*x13;
				T x17 = 2*f[0];
				T x18 = f[0]*f[1];
				T x19 = f[2]*f[3];
				T x20 = x18 + x19;
				T x21 = 2*f[1];
				T x22 = -x17*x3 + x17*x6 - x20*x21;
				T x23 = 3*x10;
				T x24 = -pow(x20, 2) - 1.0/2.0*pow(x3, 2) - 1.0/2.0*pow(x5, 2) + (1.0/2.0)*pow(x6, 2);
				T x25 = (15.0/4.0)*x24/pow(x9, 7.0/2.0);
				T x26 = f[2]*x13;
				T x27 = f[0]*x26;
				T x28 = f[3]*x13;
				T x29 = f[1]*x28;
				T x30 = x13*x19;
				T x31 = (3.0/2.0)*x10;
				T x32 = x22*x31;
				T x33 = -x17*x20 - x21*x5 + x21*x6;
				T x34 = f[3]*x31;
				T x35 = C1*(-x11*x19 + x27 - x29) + C2*(f[2]*x32 - x19*x25 - 2*x30 - x33*x34);
				T x36 = 2*x8;
				T x37 = x12 + x13*x36;
				T x38 = x13*x2;
				T x39 = x13*x18;
				T x40 = f[1]*f[3];
				T x41 = 2*f[2];
				T x42 = 2*f[3];
				T x43 = -x20*x42 - x3*x41 + x41*x6;
				T x44 = grad_test[0]*(C1*(-x11*x40 - x30 + x39) + C2*(f[1]*x32 - x21*x28 - x25*x40 - x34*x43));
				T x45 = x13*x4;
				T x46 = (1.0/2.0)*x13*x6;
				T x47 = x31*x33;
				T x48 = f[2]*x31;
				T x49 = x24*x31;
				T x50 = C1*(x11*x8 + x38 + x45 + x46) + C2*(f[1]*x47 + x13*(-x14 + 4*x8) + x25*x8 + x43*x48 + x49);
				T x51 = f[0]*f[2];
				T x52 = -x20*x41 - x42*x5 + x42*x6;
				T x53 = grad_test[1]*(C1*(-x11*x51 + x30 - x39) + C2*(-f[0]*x47 - x17*x26 - x25*x51 + x48*x52));
				T x54 = x1*x13;
				T x55 = C1*(x11*x7 - x16 - x46 - x54) + C2*(-f[0]*x32 + x13*(-x36 + 4*x7) + x25*x7 - x34*x52 - x49);
				T x56 = C1*(-x11*x18 - x27 + x29) + C2*(-f[0]*x31*x43 + f[1]*x31*x52 - x18*x25 - 2*x39);
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(C1*(x0*x11 + x15) + C2*(-f[3]*x22*x23 + x0*x25 + 2*x16)) + grad_test[1]*x35) + grad_trial[1]*(grad_test[0]*x35 + grad_test[1]*(C1*(x11*x2 + x37) + C2*(f[2]*x23*x33 + x2*x25 + 2*x38))));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x50 + x44) + grad_trial[1]*(grad_test[0]*x55 + x53));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x55 + x44) + grad_trial[1]*(grad_test[0]*x50 + x53));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(C1*(x11*x4 + x37) + C2*(f[1]*x23*x43 + x25*x4 + 2*x45)) + grad_test[1]*x56) + grad_trial[1]*(grad_test[0]*x56 + grad_test[1]*(C1*(x1*x11 + x15) + C2*(-f[0]*x23*x52 + x1*x25 + 2*x54))));
			}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T *UTOPIA_RESTRICT grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = f[0]*f[3] - f[1]*f[2];
				T x1 = pow(x0, -1.0/2.0);
				T x2 = 2*f[0];
				T x3 = pow(x0, -3.0/2.0);
				T x4 = pow(f[0], 2) + pow(f[2], 2);
				T x5 = pow(f[1], 2) + pow(f[3], 2);
				T x6 = x4 + x5;
				T x7 = (1.0/2.0)*x3*x6;
				T x8 = f[0]*f[1] + f[2]*f[3];
				T x9 = 2*f[1];
				T x10 = (3.0/2.0)*(-1.0/2.0*pow(x4, 2) - 1.0/2.0*pow(x5, 2) + (1.0/2.0)*pow(x6, 2) - pow(x8, 2))/pow(x0, 5.0/2.0);
				T x11 = 2*f[2];
				T x12 = 2*f[3];
				lf[0] += dx*(grad_test[0]*(C1*(-f[3]*x7 + x1*x2) + C2*(-f[3]*x10 + x3*(-x2*x4 + x2*x6 - x8*x9))) + grad_test[1]*(C1*(f[2]*x7 + x1*x9) + C2*(f[2]*x10 + x3*(-x2*x8 - x5*x9 + x6*x9))));
				lf[1] += dx*(grad_test[0]*(C1*(f[1]*x7 + x1*x11) + C2*(f[1]*x10 + x3*(-x11*x4 + x11*x6 - x12*x8))) + grad_test[1]*(C1*(-f[0]*x7 + x1*x12) + C2*(-f[0]*x10 + x3*(-x11*x8 - x12*x5 + x12*x6))));
			}

			UTOPIA_FUNCTION void value(
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = f[0]*f[3] - f[1]*f[2];
				T x1 = pow(f[0], 2) + pow(f[2], 2);
				T x2 = pow(f[1], 2) + pow(f[3], 2);
				T x3 = x1 + x2;
				e += dx*(C1*(-2 + x3/sqrt(x0)) + C2*(-2 + (-1.0/2.0*pow(x1, 2) - 1.0/2.0*pow(x2, 2) + (1.0/2.0)*pow(x3, 2) - pow(f[0]*f[1] + f[2]*f[3], 2))/pow(x0, 3.0/2.0)));
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
				T x0 = f[0]*f[3];
				T x1 = f[1]*f[2];
				T x2 = x0 - x1;
				T x3 = pow(x2, -1.0/2.0);
				T x4 = pow(f[0], 2);
				T x5 = pow(f[2], 2);
				T x6 = x4 + x5;
				T x7 = pow(f[1], 2);
				T x8 = pow(f[3], 2);
				T x9 = x7 + x8;
				T x10 = x6 + x9;
				T x11 = pow(x2, -3.0/2.0);
				T x12 = f[0]*f[1];
				T x13 = f[2]*f[3];
				T x14 = x12 + x13;
				T x15 = (1.0/2.0)*pow(x10, 2) - pow(x14, 2) - 1.0/2.0*pow(x6, 2) - 1.0/2.0*pow(x9, 2);
				T x16 = 2*x3;
				T x17 = f[3]*x11;
				T x18 = (1.0/2.0)*x10;
				T x19 = 2*f[0];
				T x20 = 2*f[1];
				T x21 = x10*x19 - x14*x20 - x19*x6;
				T x22 = pow(x2, -5.0/2.0);
				T x23 = (3.0/2.0)*x22;
				T x24 = x15*x23;
				T x25 = f[2]*x11;
				T x26 = x10*x20 - x14*x19 - x20*x9;
				T x27 = x11*x18;
				T x28 = 2*f[2];
				T x29 = 2*f[3];
				T x30 = x10*x28 - x14*x29 - x28*x6;
				T x31 = x10*x29 - x14*x28 - x29*x9;
				T x32 = (3.0/4.0)*x10*x22;
				T x33 = 2*x0;
				T x34 = -x11*x33 + x16;
				T x35 = x11*x8;
				T x36 = 3*x22;
				T x37 = (15.0/4.0)*x15/pow(x2, 7.0/2.0);
				T x38 = f[0]*x25;
				T x39 = f[1]*x17;
				T x40 = x11*x13;
				T x41 = x21*x23;
				T x42 = f[3]*x23;
				T x43 = C1*(-x13*x32 + x38 - x39) + C2*(f[2]*x41 - x13*x37 - x26*x42 - 2*x40);
				T x44 = 2*x1;
				T x45 = x11*x44 + x16;
				T x46 = x11*x5;
				T x47 = x11*x12;
				T x48 = f[1]*f[3];
				T x49 = grad_test[0]*(C1*(-x32*x48 - x40 + x47) + C2*(f[1]*x41 - x17*x20 - x30*x42 - x37*x48));
				T x50 = x11*x7;
				T x51 = x23*x26;
				T x52 = f[2]*x23;
				T x53 = C1*(x1*x32 + x27 + x46 + x50) + C2*(f[1]*x51 + x1*x37 + x11*(4*x1 - x33) + x24 + x30*x52);
				T x54 = f[0]*f[2];
				T x55 = grad_test[1]*(C1*(-x32*x54 + x40 - x47) + C2*(-f[0]*x51 - x19*x25 + x31*x52 - x37*x54));
				T x56 = x11*x4;
				T x57 = C1*(x0*x32 - x27 - x35 - x56) + C2*(-f[0]*x41 + x0*x37 + x11*(4*x0 - x44) - x24 - x31*x42);
				T x58 = C1*(-x12*x32 - x38 + x39) + C2*(-f[0]*x23*x30 + f[1]*x23*x31 - x12*x37 - 2*x47);
				e += dx*(C1*(x10*x3 - 2) + C2*(x11*x15 - 2));
				lf[0] += dx*(grad_test[0]*(C1*(f[0]*x16 - x17*x18) + C2*(-f[3]*x24 + x11*x21)) + grad_test[1]*(C1*(f[1]*x16 + x18*x25) + C2*(f[2]*x24 + x11*x26)));
				lf[1] += dx*(grad_test[0]*(C1*(f[1]*x27 + f[2]*x16) + C2*(f[1]*x24 + x11*x30)) + grad_test[1]*(C1*(-f[0]*x27 + f[3]*x16) + C2*(-f[0]*x24 + x11*x31)));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(C1*(x32*x8 + x34) + C2*(-f[3]*x21*x36 + 2*x35 + x37*x8)) + grad_test[1]*x43) + grad_trial[1]*(grad_test[0]*x43 + grad_test[1]*(C1*(x32*x5 + x45) + C2*(f[2]*x26*x36 + x37*x5 + 2*x46))));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x53 + x49) + grad_trial[1]*(grad_test[0]*x57 + x55));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x57 + x49) + grad_trial[1]*(grad_test[0]*x53 + x55));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(C1*(x32*x7 + x45) + C2*(f[1]*x30*x36 + x37*x7 + 2*x50)) + grad_test[1]*x58) + grad_trial[1]*(grad_test[0]*x58 + grad_test[1]*(C1*(x32*x4 + x34) + C2*(-f[0]*x31*x36 + x37*x4 + 2*x56))));
			}

			T C1{1.0};
			T C2{1.0};

		};
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_2_IMPL_hpp
