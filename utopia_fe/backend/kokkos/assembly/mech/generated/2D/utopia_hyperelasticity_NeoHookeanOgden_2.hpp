#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanOgden.hpp"

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of NeoHookeanOgden for dimension 2
		 */
		template<typename T>
		class NeoHookeanOgden<T, 2> {
		public:
			static constexpr int Dim = 2;

			class Params : public Configurable {
			public:
				void read(Input &in) override
				{
					in.get("mu", mu);
					in.get("lambda", lambda);
				}

				T mu{1.0};
				T lambda{1.0};
			};

			NeoHookeanOgden(const Params &params)
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
				T x0 = f[0]*f[3];
				T x1 = f[1]*f[2];
				T x2 = x0 - x1;
				T x3 = pow(x2, -2);
				T x4 = f[3]*x3;
				T x5 = f[2]*x4;
				T x6 = log(x2);
				T x7 = lambda*x6;
				T x8 = -lambda*x5 - mu*x5 + x5*x7;
				T x9 = pow(f[3], 2)*x3;
				T x10 = lambda*x9;
				T x11 = pow(f[2], 2)*x3;
				T x12 = lambda*x11;
				T x13 = f[1]*x4;
				T x14 = grad_test[0]*(-lambda*x13 - mu*x13 + x13*x7);
				T x15 = 1.0/x2;
				T x16 = mu*x15;
				T x17 = x1*x3;
				T x18 = x15*x7;
				T x19 = lambda*x17 + mu*x17 + x16 - x17*x7 - x18;
				T x20 = f[0]*x3;
				T x21 = f[2]*x20;
				T x22 = grad_test[1]*(-lambda*x21 - mu*x21 + x21*x7);
				T x23 = x0*x3;
				T x24 = lambda*x23 + mu*x23 - x16 + x18 - x23*x7;
				T x25 = f[1]*x20;
				T x26 = -lambda*x25 - mu*x25 + x25*x7;
				T x27 = pow(f[1], 2)*x3;
				T x28 = pow(f[0], 2)*x3;
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(mu*x9 + mu - x10*x6 + x10) + grad_test[1]*x8) + grad_trial[1]*(grad_test[0]*x8 + grad_test[1]*(mu*x11 + mu - x12*x6 + x12)));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x19 + x14) + grad_trial[1]*(grad_test[0]*x24 + x22));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x24 + x14) + grad_trial[1]*(grad_test[0]*x19 + x22));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(lambda*x27 + mu*x27 + mu - x27*x7) + grad_test[1]*x26) + grad_trial[1]*(grad_test[0]*x26 + grad_test[1]*(lambda*x28 + mu*x28 + mu - x28*x7)));
			}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T *UTOPIA_RESTRICT grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = f[0]*mu;
				T x1 = f[0]*f[3] - f[1]*f[2];
				T x2 = 1.0/x1;
				T x3 = f[3]*mu;
				T x4 = lambda*x2*log(x1);
				T x5 = f[1]*mu;
				T x6 = f[2]*mu;
				lf[0] += dx*(grad_test[0]*(f[3]*x4 + x0 - x2*x3) + grad_test[1]*(-f[2]*x4 + x2*x6 + x5));
				lf[1] += dx*(grad_test[0]*(-f[1]*x4 + x2*x5 + x6) + grad_test[1]*(f[0]*x4 - x0*x2 + x3));
			}

			UTOPIA_FUNCTION void value(
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = log(f[0]*f[3] - f[1]*f[2]);
				e += dx*((1.0/2.0)*lambda*pow(x0, 2) - mu*x0 + (1.0/2.0)*mu*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) - 2));
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
				T x3 = log(x2);
				T x4 = pow(f[0], 2);
				T x5 = pow(f[1], 2);
				T x6 = pow(f[2], 2);
				T x7 = pow(f[3], 2);
				T x8 = f[0]*mu;
				T x9 = 1.0/x2;
				T x10 = f[3]*mu;
				T x11 = lambda*x3;
				T x12 = x11*x9;
				T x13 = f[1]*mu;
				T x14 = f[2]*mu;
				T x15 = pow(x2, -2);
				T x16 = f[3]*x15;
				T x17 = lambda*x16;
				T x18 = x11*x16;
				T x19 = -f[2]*x17 + f[2]*x18 - x14*x16;
				T x20 = x15*x7;
				T x21 = lambda*x20;
				T x22 = x15*x6;
				T x23 = lambda*x22;
				T x24 = grad_test[0]*(-f[1]*x17 + f[1]*x18 - x13*x16);
				T x25 = mu*x9;
				T x26 = x1*x15;
				T x27 = lambda*x26 + mu*x26 - x11*x26 - x12 + x25;
				T x28 = f[2]*x15;
				T x29 = f[0]*lambda;
				T x30 = f[0]*x11;
				T x31 = grad_test[1]*(-x28*x29 + x28*x30 - x28*x8);
				T x32 = x0*x15;
				T x33 = lambda*x32 + mu*x32 - x11*x32 + x12 - x25;
				T x34 = f[1]*x15;
				T x35 = -x29*x34 + x30*x34 - x34*x8;
				T x36 = x15*x5;
				T x37 = x15*x4;
				e += dx*((1.0/2.0)*lambda*pow(x3, 2) - mu*x3 + (1.0/2.0)*mu*(x4 + x5 + x6 + x7 - 2));
				lf[0] += dx*(grad_test[0]*(f[3]*x12 - x10*x9 + x8) + grad_test[1]*(-f[2]*x12 + x13 + x14*x9));
				lf[1] += dx*(grad_test[0]*(-f[1]*x12 + x13*x9 + x14) + grad_test[1]*(f[0]*x12 + x10 - x8*x9));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(mu*x20 + mu - x21*x3 + x21) + grad_test[1]*x19) + grad_trial[1]*(grad_test[0]*x19 + grad_test[1]*(mu*x22 + mu - x23*x3 + x23)));
				bf[1] += dx*(grad_trial[0]*(grad_test[1]*x27 + x24) + grad_trial[1]*(grad_test[0]*x33 + x31));
				bf[2] += dx*(grad_trial[0]*(grad_test[1]*x33 + x24) + grad_trial[1]*(grad_test[0]*x27 + x31));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(lambda*x36 + mu*x36 + mu - x11*x36) + grad_test[1]*x35) + grad_trial[1]*(grad_test[0]*x35 + grad_test[1]*(lambda*x37 + mu*x37 + mu - x11*x37)));
			}

			T mu{1.0};
			T lambda{1.0};

		};
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_2_IMPL_hpp
