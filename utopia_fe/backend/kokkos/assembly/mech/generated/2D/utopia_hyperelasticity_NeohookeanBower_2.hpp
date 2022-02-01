#ifndef UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_2_IMPL_hpp

#include "utopia_hyperelasticity_NeohookeanBower.hpp"

namespace utopia {
	namespace kernels {

		/** 
		 * Specialization of NeohookeanBower for dimension 2 
		 */
		template<typename T>
		class NeohookeanBower<T, 2> {
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
				T x14 = f[3]*lmbda;
				T x15 = (4.0/3.0)*x10;
				T x16 = f[0]*x15;
				T x17 = f[2]*x16;
				T x18 = f[3]*x15;
				T x19 = f[1]*x18;
				T x20 = f[3]*x8;
				T x21 = -f[2]*x14 + x13*(-f[2]*x20 + x17 - x19);
				T x22 = x11*x6 + x9;
				T x23 = f[0]*lmbda;
				T x24 = f[1]*x16;
				T x25 = f[2]*x18;
				T x26 = f[0]*x8;
				T x27 = grad_test[1]*(-f[2]*x23 + x13*(-f[2]*x26 - x24 + x25));
				T x28 = lmbda*(x7 - 1);
				T x29 = (2.0/3.0)*x10*x4;
				T x30 = lmbda*x5 + x13*(-x0*x15 - x1*x15 - x29 + x5*x8) + x28;
				T x31 = grad_test[0]*(-f[1]*x14 + x13*(-f[1]*x20 + x24 - x25));
				T x32 = lmbda*x6 + x13*(x15*x2 + x15*x3 + x29 + x6*x8) - x28;
				T x33 = -f[1]*x23 + x13*(-f[1]*x26 - x17 + x19);
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(lmbda*x0 + x13*(x0*x8 + x12)) + grad_test[1]*x21) + grad_trial[1]*(grad_test[0]*x21 + grad_test[1]*(lmbda*x3 + x13*(x22 + x3*x8))));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x32 + x31) + grad_trial[1]*(grad_test[0]*x30 + x27));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x30 + x31) + grad_trial[1]*(grad_test[0]*x32 + x27));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(lmbda*x2 + x13*(x2*x8 + x22)) + grad_test[1]*x33) + grad_trial[1]*(grad_test[0]*x33 + grad_test[1]*(lmbda*x1 + x13*(x1*x8 + x12))));
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
				T x0 = f[0]*f[3] - f[1]*f[2];
				T x1 = lmbda*(x0 - 1);
				T x2 = 2/pow(x0, 2.0/3.0);
				T x3 = (2.0/3.0)*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2))/pow(x0, 5.0/3.0);
				T x4 = (1.0/2.0)*mu;
				lf[0] += dx*(grad_test[0]*(f[3]*x1 + x4*(f[0]*x2 - f[3]*x3)) + grad_test[1]*(-f[2]*x1 + x4*(f[1]*x2 + f[2]*x3)));
				lf[1] += dx*(grad_test[0]*(-f[1]*x1 + x4*(f[1]*x3 + f[2]*x2)) + grad_test[1]*(f[0]*x1 + x4*(-f[0]*x3 + f[3]*x2)));
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
				T x0 = f[0]*f[3] - f[1]*f[2];
				e += dx*((1.0/2.0)*lmbda*pow(x0 - 1, 2) + (1.0/2.0)*mu*(-2 + (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2))/pow(x0, 2.0/3.0)));
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
				T x0 = f[0]*f[3];
				T x1 = f[1]*f[2];
				T x2 = x0 - x1;
				T x3 = x2 - 1;
				T x4 = pow(x2, -2.0/3.0);
				T x5 = pow(f[0], 2);
				T x6 = pow(f[1], 2);
				T x7 = pow(f[2], 2);
				T x8 = pow(f[3], 2);
				T x9 = x5 + x6 + x7 + x8;
				T x10 = (1.0/2.0)*mu;
				T x11 = lmbda*x3;
				T x12 = 2*x4;
				T x13 = pow(x2, -5.0/3.0);
				T x14 = (2.0/3.0)*x13*x9;
				T x15 = (10.0/9.0)*x9/pow(x2, 8.0/3.0);
				T x16 = (8.0/3.0)*x13;
				T x17 = -x0*x16 + x12;
				T x18 = f[3]*lmbda;
				T x19 = (4.0/3.0)*x13;
				T x20 = f[0]*x19;
				T x21 = f[2]*x20;
				T x22 = f[3]*x19;
				T x23 = f[1]*x22;
				T x24 = f[3]*x15;
				T x25 = -f[2]*x18 + x10*(-f[2]*x24 + x21 - x23);
				T x26 = x1*x16 + x12;
				T x27 = f[0]*lmbda;
				T x28 = f[1]*x20;
				T x29 = f[2]*x22;
				T x30 = f[0]*x15;
				T x31 = grad_test[1]*(-f[2]*x27 + x10*(-f[2]*x30 - x28 + x29));
				T x32 = lmbda*x0 + x10*(x0*x15 - x14 - x19*x5 - x19*x8) + x11;
				T x33 = grad_test[0]*(-f[1]*x18 + x10*(-f[1]*x24 + x28 - x29));
				T x34 = lmbda*x1 + x10*(x1*x15 + x14 + x19*x6 + x19*x7) - x11;
				T x35 = -f[1]*x27 + x10*(-f[1]*x30 - x21 + x23);
				e += dx*((1.0/2.0)*lmbda*pow(x3, 2) + x10*(x4*x9 - 2));
				lf[0] += dx*(grad_test[0]*(f[3]*x11 + x10*(f[0]*x12 - f[3]*x14)) + grad_test[1]*(-f[2]*x11 + x10*(f[1]*x12 + f[2]*x14)));
				lf[1] += dx*(grad_test[0]*(-f[1]*x11 + x10*(f[1]*x14 + f[2]*x12)) + grad_test[1]*(f[0]*x11 + x10*(-f[0]*x14 + f[3]*x12)));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(lmbda*x8 + x10*(x15*x8 + x17)) + grad_test[1]*x25) + grad_trial[1]*(grad_test[0]*x25 + grad_test[1]*(lmbda*x7 + x10*(x15*x7 + x26))));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x34 + x33) + grad_trial[1]*(grad_test[0]*x32 + x31));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x32 + x33) + grad_trial[1]*(grad_test[0]*x34 + x31));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(lmbda*x6 + x10*(x15*x6 + x26)) + grad_test[1]*x35) + grad_trial[1]*(grad_test[0]*x35 + grad_test[1]*(lmbda*x5 + x10*(x15*x5 + x17))));
			}

		};
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_2_IMPL_hpp
