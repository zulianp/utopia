#ifndef UTOPIA_TPL_HYPERELASTICITY_NeohookeanSmith_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeohookeanSmith_2_IMPL_hpp

#include "utopia_hyperelasticity_NeohookeanSmith.hpp"

namespace utopia {
	namespace kernels {

		/** 
		 * Specialization of NeohookeanSmith for dimension 2 
		 */
		template<typename T>
		class NeohookeanSmith<T, 2> {
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
				T x0 = (11.0/6.0)*lmbda;
				T x1 = f[3]*x0;
				T x2 = pow(f[0], 2);
				T x3 = pow(f[1], 2);
				T x4 = pow(f[2], 2);
				T x5 = pow(f[3], 2);
				T x6 = x2 + x3 + x4 + x5 + 1;
				T x7 = (8.0/3.0)*mu/pow(x6, 2);
				T x8 = f[0]*x7;
				T x9 = f[1]*x8 - f[2]*x1;
				T x10 = (4.0/3.0)*mu;
				T x11 = x10 - x10/x6;
				T x12 = grad_test[0]*(-f[1]*x1 + f[2]*x8);
				T x13 = f[1]*f[2];
				T x14 = f[0]*f[3];
				T x15 = x0*(-x13 + x14 - 1 - 6.0/11.0*mu/lmbda);
				T x16 = x0*x13 + x13*x7 - x15;
				T x17 = f[0]*x0;
				T x18 = f[3]*x7;
				T x19 = grad_test[1]*(f[1]*x18 - f[2]*x17);
				T x20 = x0*x14 + x14*x7 + x15;
				T x21 = -f[1]*x17 + f[2]*x18;
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x0*x5 + x11 + x2*x7) + grad_test[1]*x9) + grad_trial[1]*(grad_test[0]*x9 + grad_test[1]*(x0*x4 + x11 + x3*x7)));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x16 + x12) + grad_trial[1]*(grad_test[0]*x20 + x19));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x20 + x12) + grad_trial[1]*(grad_test[0]*x16 + x19));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x0*x3 + x11 + x4*x7) + grad_test[1]*x21) + grad_trial[1]*(grad_test[0]*x21 + grad_test[1]*(x0*x2 + x11 + x5*x7)));
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
				T x0 = (4.0/3.0)*mu;
				T x1 = f[0]*x0;
				T x2 = (11.0/6.0)*lmbda*(f[0]*f[3] - f[1]*f[2] - 1 - 6.0/11.0*mu/lmbda);
				T x3 = 1.0/(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + 1);
				T x4 = f[1]*x0;
				T x5 = f[2]*x0;
				T x6 = f[3]*x0;
				lf[0] += dx*(grad_test[0]*(f[3]*x2 - x1*x3 + x1) + grad_test[1]*(-f[2]*x2 - x3*x4 + x4));
				lf[1] += dx*(grad_test[0]*(-f[1]*x2 - x3*x5 + x5) + grad_test[1]*(f[0]*x2 - x3*x6 + x6));
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
				T x0 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2);
				T x1 = (2.0/3.0)*mu;
				e += dx*((11.0/12.0)*lmbda*pow(f[0]*f[3] - f[1]*f[2] - 1 - 6.0/11.0*mu/lmbda, 2) + x1*(x0 - 2) - x1*log(x0 + 1));
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
				T x0 = pow(f[0], 2);
				T x1 = pow(f[1], 2);
				T x2 = pow(f[2], 2);
				T x3 = pow(f[3], 2);
				T x4 = x0 + x1 + x2 + x3;
				T x5 = (2.0/3.0)*mu;
				T x6 = x4 + 1;
				T x7 = f[0]*f[3];
				T x8 = f[1]*f[2];
				T x9 = x7 - x8 - 1 - 6.0/11.0*mu/lmbda;
				T x10 = (4.0/3.0)*mu;
				T x11 = f[0]*x10;
				T x12 = (11.0/6.0)*lmbda;
				T x13 = x12*x9;
				T x14 = 1.0/x6;
				T x15 = f[1]*x10;
				T x16 = f[2]*x10;
				T x17 = f[3]*x10;
				T x18 = f[3]*x12;
				T x19 = (8.0/3.0)*mu/pow(x6, 2);
				T x20 = f[0]*x19;
				T x21 = f[1]*x20 - f[2]*x18;
				T x22 = -x10*x14 + x10;
				T x23 = grad_test[0]*(-f[1]*x18 + f[2]*x20);
				T x24 = x12*x8 - x13 + x19*x8;
				T x25 = f[0]*x12;
				T x26 = f[3]*x19;
				T x27 = grad_test[1]*(f[1]*x26 - f[2]*x25);
				T x28 = x12*x7 + x13 + x19*x7;
				T x29 = -f[1]*x25 + f[2]*x26;
				e += dx*((11.0/12.0)*lmbda*pow(x9, 2) + x5*(x4 - 2) - x5*log(x6));
				lf[0] += dx*(grad_test[0]*(f[3]*x13 - x11*x14 + x11) + grad_test[1]*(-f[2]*x13 - x14*x15 + x15));
				lf[1] += dx*(grad_test[0]*(-f[1]*x13 - x14*x16 + x16) + grad_test[1]*(f[0]*x13 - x14*x17 + x17));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x0*x19 + x12*x3 + x22) + grad_test[1]*x21) + grad_trial[1]*(grad_test[0]*x21 + grad_test[1]*(x1*x19 + x12*x2 + x22)));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x24 + x23) + grad_trial[1]*(grad_test[0]*x28 + x27));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x28 + x23) + grad_trial[1]*(grad_test[0]*x24 + x27));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x1*x12 + x19*x2 + x22) + grad_test[1]*x29) + grad_trial[1]*(grad_test[0]*x29 + grad_test[1]*(x0*x12 + x19*x3 + x22)));
			}

		};;
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_NeohookeanSmith_2_IMPL_hpp
