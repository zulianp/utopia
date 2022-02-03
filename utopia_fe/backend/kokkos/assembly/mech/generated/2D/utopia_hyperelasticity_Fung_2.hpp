#ifndef UTOPIA_TPL_HYPERELASTICITY_Fung_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_Fung_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_Fung.hpp"

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of Fung for dimension 2
		 */
		template<typename T>
		class Fung<T, 2> {
		public:
			static constexpr int Dim = 2;

			class Params : public Configurable {
			public:
				void read(Input &in) override
				{
					in.get("a", a);
					in.get("b", b);
					in.get("c", c);
				}

				T a{1.0};
				T b{1.0};
				T c{1.0};
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
				T x4 = b*exp(c*(x0 + x1 + x2 + x3 - 2) - 1);
				T x5 = 2*pow(c, 2)*x4;
				T x6 = f[0]*x5;
				T x7 = f[1]*x6;
				T x8 = a + c*x4;
				T x9 = grad_test[0]*x6;
				T x10 = f[2]*x9;
				T x11 = f[2]*x5;
				T x12 = f[1]*grad_test[1];
				T x13 = f[3]*x12*x5;
				T x14 = f[3]*grad_test[1];
				T x15 = grad_test[0]*x11;
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x0*x5 + x8) + grad_test[1]*x7) + grad_trial[1]*(grad_test[0]*x7 + grad_test[1]*(x1*x5 + x8)));
				bf[1] += dx*(grad_trial[0]*(x10 + x11*x12) + grad_trial[1]*(f[3]*x9 + x13));
				bf[2] += dx*(grad_trial[0]*(x10 + x14*x6) + grad_trial[1]*(f[1]*x15 + x13));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x2*x5 + x8) + x11*x14) + grad_trial[1]*(f[3]*x15 + grad_test[1]*(x3*x5 + x8)));
			}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T *UTOPIA_RESTRICT grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = b*c*exp(c*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) - 2) - 1);
				lf[0] += dx*(grad_test[0]*(a*f[0] + f[0]*x0) + grad_test[1]*(a*f[1] + f[1]*x0));
				lf[1] += dx*(grad_test[0]*(a*f[2] + f[2]*x0) + grad_test[1]*(a*f[3] + f[3]*x0));
			}

			UTOPIA_FUNCTION void value(
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) - 2;
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
				T x4 = x0 + x1 + x2 + x3 - 2;
				T x5 = b*exp(c*x4 - 1);
				T x6 = c*x5;
				T x7 = 2*pow(c, 2)*x5;
				T x8 = f[0]*x7;
				T x9 = f[1]*x8;
				T x10 = a + x6;
				T x11 = grad_test[0]*x8;
				T x12 = f[2]*x11;
				T x13 = f[2]*x7;
				T x14 = f[1]*grad_test[1];
				T x15 = f[3]*x14*x7;
				T x16 = f[3]*grad_test[1];
				T x17 = grad_test[0]*x13;
				e += dx*((1.0/2.0)*a*x4 + (1.0/2.0)*x5);
				lf[0] += dx*(grad_test[0]*(a*f[0] + f[0]*x6) + grad_test[1]*(a*f[1] + f[1]*x6));
				lf[1] += dx*(grad_test[0]*(a*f[2] + f[2]*x6) + grad_test[1]*(a*f[3] + f[3]*x6));
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x0*x7 + x10) + grad_test[1]*x9) + grad_trial[1]*(grad_test[0]*x9 + grad_test[1]*(x1*x7 + x10)));
				bf[1] += dx*(grad_trial[0]*(x12 + x13*x14) + grad_trial[1]*(f[3]*x11 + x15));
				bf[2] += dx*(grad_trial[0]*(x12 + x16*x8) + grad_trial[1]*(f[1]*x17 + x15));
				bf[3] += dx*(grad_trial[0]*(grad_test[0]*(x10 + x2*x7) + x13*x16) + grad_trial[1]*(f[3]*x17 + grad_test[1]*(x10 + x3*x7)));
			}

			T a{1.0};
			T b{1.0};
			T c{1.0};

		};
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_Fung_2_IMPL_hpp
