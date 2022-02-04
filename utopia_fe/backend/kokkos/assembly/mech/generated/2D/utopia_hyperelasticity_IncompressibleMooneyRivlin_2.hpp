#ifndef UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_IncompressibleMooneyRivlin.hpp"

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of IncompressibleMooneyRivlin for dimension 2
		 */
		template<typename T>
		class IncompressibleMooneyRivlin<T, 2> {
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

			IncompressibleMooneyRivlin(const Params &params)
			{
				C1 = params.C1;
				C2 = params.C2;
			}

			UTOPIA_FUNCTION void hessian(
				const T *UTOPIA_RESTRICT f,
				const T p,
				const T *grad_test,
				const T *grad_trial,
				const T *fun_test,
				const T *fun_trial,
				const T dx,
				T *UTOPIA_RESTRICT bf) const
			{
				using namespace utopia::device;
				// Automatically generated
				T x0 = 2*C2;
				T x1 = p + x0;
				T x2 = f[2]*grad_test[1];
				T x3 = x1*x2;
				T x4 = 2*C1;
				T x5 = pow(f[3], 2);
				T x6 = f[3]*grad_test[0];
				T x7 = x1*x6;
				T x8 = pow(f[2], 2);
				T x9 = f[1]*x7;
				T x10 = f[1]*f[2];
				T x11 = f[0]*f[3];
				T x12 = x10 - x11 + 1;
				T x13 = p*x12;
				T x14 = p*x10 - x0*(-2*x10 + x11) + x13;
				T x15 = -f[0]*x3;
				T x16 = p*x11 + x0*(-x10 + 2*x11) - x13;
				T x17 = f[0]*grad_test[1];
				T x18 = pow(f[1], 2);
				T x19 = f[1]*grad_test[0];
				T x20 = pow(f[0], 2);
				T x21 = dx*x12;
				T x22 = fun_test*x21;
				T x23 = fun_trial*x21;
				bf[0] += dx*(grad_trial[0]*(-f[3]*x3 + grad_test[0]*(p*x5 + x0*x5 + x4)) - grad_trial[1]*(f[2]*x7 - grad_test[1]*(p*x8 + x0*x8 + x4)));
				bf[1] += dx*(-grad_trial[0]*(-grad_test[1]*x14 + x9) + grad_trial[1]*(grad_test[0]*x16 + x15));
				bf[3] += dx*(-grad_trial[0]*(-grad_test[1]*x16 + x9) + grad_trial[1]*(grad_test[0]*x14 + x15));
				bf[4] += dx*(grad_trial[0]*(-f[1]*x1*x17 + grad_test[0]*(p*x18 + x0*x18 + x4)) - grad_trial[1]*(f[0]*x1*x19 - grad_test[1]*(p*x20 + x0*x20 + x4)));
				bf[8] += 0;
				bf[6] += x22*(f[2]*grad_trial[1] - f[3]*grad_trial[0]);
				bf[7] += -x22*(f[0]*grad_trial[1] - f[1]*grad_trial[0]);
				bf[2] += x23*(x2 - x6);
				bf[5] += -x23*(x17 - x19);
			}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T p,
				const T *UTOPIA_RESTRICT grad_test,
				const T *fun_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = grad_test[1]*p;
				T x1 = f[2]*x0;
				T x2 = 2*C1;
				T x3 = grad_test[0]*x2;
				T x4 = grad_test[1]*x2;
				T x5 = grad_test[0]*p;
				T x6 = f[3]*x5;
				T x7 = f[0]*pow(f[3], 2);
				T x8 = f[1]*pow(f[2], 2);
				T x9 = f[0]*f[3];
				T x10 = 2*C2;
				T x11 = grad_test[1]*x10;
				T x12 = f[1]*f[2];
				T x13 = grad_test[0]*x10;
				T x14 = f[1]*x5;
				T x15 = f[0]*x0;
				T x16 = pow(f[0], 2)*f[3];
				T x17 = pow(f[1], 2)*f[2];
				lf[0] += dx*(f[0]*x3 + f[1]*x4 - f[2]*x11*x9 - f[3]*x12*x13 + x0*x8 - x1*x9 + x1 + x11*x8 - x12*x6 + x13*x7 + x5*x7 - x6);
				lf[1] += dx*(-f[0]*x11*x12 - f[1]*x13*x9 + f[2]*x3 + f[3]*x4 + x0*x16 + x11*x16 - x12*x15 + x13*x17 - x14*x9 + x14 - x15 + x17*x5);
				lf[2] += (1.0/2.0)*dx*fun_test*pow(x12 - x9 + 1, 2);
			}

			UTOPIA_FUNCTION void value(
				const T *UTOPIA_RESTRICT f,
				const T p,
				const T dx,
				T &e
				) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = pow(f[0], 2) + pow(f[2], 2);
				T x1 = pow(f[1], 2) + pow(f[3], 2);
				T x2 = x0 + x1;
				e += dx*(C1*(x2 - 2) + C2*(-1.0/2.0*pow(x0, 2) - 1.0/2.0*pow(x1, 2) + (1.0/2.0)*pow(x2, 2) - pow(f[0]*f[1] + f[2]*f[3], 2) - 2) + (1.0/2.0)*p*pow(f[0]*f[3] - f[1]*f[2] - 1, 2));
			}

			UTOPIA_FUNCTION void eval(
				const T *UTOPIA_RESTRICT f,
				const T p,
				const T *grad_test,
				const T *grad_trial,
				const T *fun_test,
				const T *fun_trial,
				const T dx,
				T &e,
				T *UTOPIA_RESTRICT lf,
				T *UTOPIA_RESTRICT bf) const
			{
				using namespace utopia::device;
			    // Automatically generated
				T x0 = f[0]*f[3];
				T x1 = f[1]*f[2];
				T x2 = -x1;
				T x3 = pow(f[0], 2);
				T x4 = pow(f[2], 2);
				T x5 = x3 + x4;
				T x6 = pow(f[1], 2);
				T x7 = pow(f[3], 2);
				T x8 = x6 + x7;
				T x9 = x5 + x8;
				T x10 = f[0]*f[1];
				T x11 = f[2]*f[3];
				T x12 = f[2]*grad_test[1];
				T x13 = p*x12;
				T x14 = 2*C1;
				T x15 = f[0]*grad_test[0];
				T x16 = f[1]*grad_test[1];
				T x17 = f[3]*grad_test[0];
				T x18 = p*x17;
				T x19 = p*x7;
				T x20 = p*x4;
				T x21 = f[0]*grad_test[1];
				T x22 = 2*C2;
				T x23 = x11*x22;
				T x24 = f[1]*grad_test[0];
				T x25 = p*x21;
				T x26 = p*x24;
				T x27 = x22*x7;
				T x28 = x22*x4;
				T x29 = f[2]*grad_test[0];
				T x30 = f[3]*grad_test[1];
				T x31 = p*x3;
				T x32 = p*x6;
				T x33 = x10*x22;
				T x34 = x22*x3;
				T x35 = x22*x6;
				T x36 = -x0 + x1 + 1;
				T x37 = dx*fun_test;
				T x38 = p + x22;
				T x39 = x11*x38;
				T x40 = f[1]*x17*x38;
				T x41 = p*x36;
				T x42 = p*x1 - x22*(x0 - 2*x1) + x41;
				T x43 = -f[0]*x12*x38;
				T x44 = p*x0 + x22*(2*x0 + x2) - x41;
				T x45 = x10*x38;
				T x46 = x36*x37;
				T x47 = dx*fun_trial*x36;
				e += dx*(C1*(x9 - 2) + C2*(-1.0/2.0*pow(x5, 2) - 1.0/2.0*pow(x8, 2) + (1.0/2.0)*pow(x9, 2) - pow(x10 + x11, 2) - 2) + (1.0/2.0)*p*pow(x0 + x2 - 1, 2));
				lf[0] += dx*(-x11*x25 - x11*x26 + x13 + x14*x15 + x14*x16 + x15*x19 + x15*x27 + x16*x20 + x16*x28 - x18 - x21*x23 - x23*x24);
				lf[1] += dx*(-x10*x13 - x10*x18 - x12*x33 + x14*x29 + x14*x30 - x17*x33 - x25 + x26 + x29*x32 + x29*x35 + x30*x31 + x30*x34);
				lf[2] += (1.0/2.0)*pow(x36, 2)*x37;
				bf[0] += dx*(grad_trial[0]*(grad_test[0]*(x14 + x19 + x27) - grad_test[1]*x39) - grad_trial[1]*(grad_test[0]*x39 - grad_test[1]*(x14 + x20 + x28)));
				bf[1] += dx*(-grad_trial[0]*(-grad_test[1]*x42 + x40) + grad_trial[1]*(grad_test[0]*x44 + x43));
				bf[3] += dx*(-grad_trial[0]*(-grad_test[1]*x44 + x40) + grad_trial[1]*(grad_test[0]*x42 + x43));
				bf[4] += dx*(grad_trial[0]*(grad_test[0]*(x14 + x32 + x35) - grad_test[1]*x45) - grad_trial[1]*(grad_test[0]*x45 - grad_test[1]*(x14 + x31 + x34)));
				bf[8] += 0;
				bf[6] += x46*(f[2]*grad_trial[1] - f[3]*grad_trial[0]);
				bf[7] += -x46*(f[0]*grad_trial[1] - f[1]*grad_trial[0]);
				bf[2] += x47*(x12 - x17);
				bf[5] += -x47*(x21 - x24);
			}

			T C1{1.0};
			T C2{1.0};

		};
	}
}

#endif // UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_2_IMPL_hpp
