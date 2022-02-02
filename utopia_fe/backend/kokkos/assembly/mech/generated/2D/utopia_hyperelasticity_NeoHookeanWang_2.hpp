#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanWang.hpp"

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of NeoHookeanWang for dimension 2
		 */
		template<typename T>
		class NeoHookeanWang<T, 2> {
		public:

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
				T x2 = pow(f[1], 2);
				T x3 = pow(f[2], 2);
				T x4 = x0 + x1 + x2 + x3;
				T x5 = f[0]*f[3];
				T x6 = f[1]*f[2];
				T x7 = x5 - x6;
				T x8 = (10.0/9.0)*x4/pow(x7, 8.0/3.0);
				T x9 = 2/pow(x7, 2.0/3.0);
				T x10 = pow(x7, -5.0/3.0);
				T x11 = (8.0/3.0)*x10;
				T x12 = -x11*x5 + x9;
				T x13 = (1.0/2.0)*mu;
				T x14 = grad_test[0]*x13;
				T x15 = (4.0/3.0)*x10;
				T x16 = f[0]*x15;
				T x17 = f[2]*x16;
				T x18 = f[3]*x15;
				T x19 = f[1]*x18;
				T x20 = f[3]*x8;
				T x21 = -f[2]*x20 + x17 - x19;
				T x22 = grad_test[1]*x13;
				T x23 = x11*x6 + x9;
				T x24 = f[1]*x16;
				T x25 = f[2]*x18;
				T x26 = x14*(-f[1]*x20 + x24 - x25);
				T x27 = (1.0/2.0)*lambda;
				T x28 = (2.0/3.0)*x10*x4;
				T x29 = x13*(x15*x2 + x15*x3 + x28 + x6*x8) - x27;
				T x30 = f[0]*x8;
				T x31 = x22*(-f[2]*x30 - x24 + x25);
				T x32 = x13*(-x0*x15 - x1*x15 - x28 + x5*x8) + x27;
				T x33 = -f[1]*x30 - x17 + x19;
				bf[0] += dx*(grad_trial[0]*(x14*(x0*x8 + x12) + x21*x22) + grad_trial[1]*(x14*x21 + x22*(x23 + x3*x8)));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x29 + x26) + grad_trial[1]*(grad_test[0]*x32 + x31));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x32 + x26) + grad_trial[1]*(grad_test[0]*x29 + x31));
				bf[3] += dx*(grad_trial[0]*(x14*(x2*x8 + x23) + x22*x33) + grad_trial[1]*(x14*x33 + x22*(x1*x8 + x12)));
			}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = (1.0/2.0)*lambda;
				T x1 = f[0]*f[3] - f[1]*f[2];
				T x2 = 2/pow(x1, 2.0/3.0);
				T x3 = (2.0/3.0)*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2))/pow(x1, 5.0/3.0);
				T x4 = (1.0/2.0)*mu;
				lf[0] += dx*(grad_test[0]*(f[3]*x0 + x4*(f[0]*x2 - f[3]*x3)) + grad_test[1]*(-f[2]*x0 + x4*(f[1]*x2 + f[2]*x3)));
				lf[1] += dx*(grad_test[0]*(-f[1]*x0 + x4*(f[1]*x3 + f[2]*x2)) + grad_test[1]*(f[0]*x0 + x4*(-f[0]*x3 + f[3]*x2)));
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
				e += dx*((1.0/2.0)*lambda*(x0 - 1) + (1.0/2.0)*mu*(-2 + (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2))/pow(x0, 2.0/3.0)));
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
				T x3 = (1.0/2.0)*lambda;
				T x4 = pow(x2, -2.0/3.0);
				T x5 = pow(f[0], 2);
				T x6 = pow(f[1], 2);
				T x7 = pow(f[2], 2);
				T x8 = pow(f[3], 2);
				T x9 = x5 + x6 + x7 + x8;
				T x10 = (1.0/2.0)*mu;
				T x11 = 2*x4;
				T x12 = pow(x2, -5.0/3.0);
				T x13 = (2.0/3.0)*x12*x9;
				T x14 = (10.0/9.0)*x9/pow(x2, 8.0/3.0);
				T x15 = (8.0/3.0)*x12;
				T x16 = -x0*x15 + x11;
				T x17 = grad_test[0]*x10;
				T x18 = (4.0/3.0)*x12;
				T x19 = f[0]*x18;
				T x20 = f[2]*x19;
				T x21 = f[3]*x18;
				T x22 = f[1]*x21;
				T x23 = f[3]*x14;
				T x24 = -f[2]*x23 + x20 - x22;
				T x25 = grad_test[1]*x10;
				T x26 = x1*x15 + x11;
				T x27 = f[1]*x19;
				T x28 = f[2]*x21;
				T x29 = x17*(-f[1]*x23 + x27 - x28);
				T x30 = x10*(x1*x14 + x13 + x18*x6 + x18*x7) - x3;
				T x31 = f[0]*x14;
				T x32 = x25*(-f[2]*x31 - x27 + x28);
				T x33 = x10*(x0*x14 - x13 - x18*x5 - x18*x8) + x3;
				T x34 = -f[1]*x31 - x20 + x22;
				e += dx*(x10*(x4*x9 - 2) + x3*(x2 - 1));
				lf[0] += dx*(grad_test[0]*(f[3]*x3 + x10*(f[0]*x11 - f[3]*x13)) + grad_test[1]*(-f[2]*x3 + x10*(f[1]*x11 + f[2]*x13)));
				lf[1] += dx*(grad_test[0]*(-f[1]*x3 + x10*(f[1]*x13 + f[2]*x11)) + grad_test[1]*(f[0]*x3 + x10*(-f[0]*x13 + f[3]*x11)));
				bf[0] += dx*(grad_trial[0]*(x17*(x14*x8 + x16) + x24*x25) + grad_trial[1]*(x17*x24 + x25*(x14*x7 + x26)));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x30 + x29) + grad_trial[1]*(grad_test[0]*x33 + x32));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x33 + x29) + grad_trial[1]*(grad_test[0]*x30 + x32));
				bf[3] += dx*(grad_trial[0]*(x17*(x14*x6 + x26) + x25*x34) + grad_trial[1]*(x17*x34 + x25*(x14*x5 + x16)));
			}

			T mu;
			T lambda;

		};

		template<typename T>
		void read_parameters(Input &in, NeoHookeanWang<T, 2> &obj)
		{
			in.get("mu", obj.mu);
			in.get("lambda", obj.lambda);
		}

	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_2_IMPL_hpp
