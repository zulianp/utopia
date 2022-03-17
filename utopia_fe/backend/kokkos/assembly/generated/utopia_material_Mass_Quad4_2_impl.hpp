#ifndef UTOPIA_TPL_MATERIAL_Mass_Quad4_2_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_Quad4_2_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Algorithms.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Quad4_2.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif


namespace utopia {
	namespace kernels {

		/**
		 * Specialization of Mass for symmetric element pair trial=test=Quad4
		 */
		template<typename T, typename GeoT>
		class Mass<Quad4<T, GeoT>> {
		public:
			using ElemT = Quad4<T, GeoT>;
			static constexpr int Dim = ElemT::Dim;

			UTOPIA_FUNCTION static constexpr const char* class_name() { return "Mass<Quad4>"; }

			class Params : public Configurable {
			public:
				void read(Input &in) override
				{
					//TODO
				}

				//TODO
			};

			Mass(const Params &params = Params())
			{
				//TODO
			}

			UTOPIA_FUNCTION void hessian(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Coefficients
				const T *UTOPIA_RESTRICT u,
				// Quadrature rule
				const T x,
				const T y,
				const T weight,
				T *UTOPIA_RESTRICT H) const
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 16*ADDAUGMENTEDASSIGNMENT + 4*MUL + 3*POW
//	- Subexpressions: 7*ADD + 46*MUL + 2*POW + 11*SUB
T x0 = 1 - x;
T x1 = 1 - y;
T x2 = weight*(-px[0]*py[1]*y + px[0]*py[1] - px[0]*py[2]*x + px[0]*py[2]*y + px[0]*py[3]*x - px[0]*py[3]*y + px[1]*py[0]*y - px[1]*py[0] + px[1]*py[2]*x - px[1]*py[3]*x + px[2]*py[0]*x - px[2]*py[0]*y - px[2]*py[1]*x - px[3]*py[0]*x + px[3]*py[0]*y + px[3]*py[1]*x);
T x3 = pow(x1, 2)*x2;
T x4 = x*x0*x3;
T x5 = x*y;
T x6 = x0*x1*x2*x5;
T x7 = 1 - x5;
T x8 = x2*x7;
T x9 = x1*x8;
T x10 = x0*x9;
T x11 = pow(x, 2);
T x12 = x11*x2;
T x13 = x1*x12*y;
T x14 = x*x9;
T x15 = x5*x8;
H[0] += pow(x0, 2)*x3;
H[1] += x4;
H[2] += x6;
H[3] += x10;
H[4] += x4;
H[5] += x11*x3;
H[6] += x13;
H[7] += x14;
H[8] += x6;
H[9] += x13;
H[10] += x12*pow(y, 2);
H[11] += x15;
H[12] += x10;
H[13] += x14;
H[14] += x15;
H[15] += x2*pow(x7, 2);
			}

			UTOPIA_FUNCTION void apply(
					// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Coefficients
				const T *UTOPIA_RESTRICT u,
				// Quadrature rule
				const T x,
				const T y,
				const T weight,
				T *UTOPIA_RESTRICT Hx) const
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 4*ADDAUGMENTEDASSIGNMENT + 4*MUL
//	- Subexpressions: 10*ADD + 39*MUL + 11*SUB
T x0 = 1 - y;
T x1 = x0*(1 - x);
T x2 = x*y;
T x3 = 1 - x2;
T x4 = x*x0;
T x5 = weight*(u[0]*x1 + u[1]*x4 + u[2]*x2 + u[3]*x3)*(-px[0]*py[1]*y + px[0]*py[1] - px[0]*py[2]*x + px[0]*py[2]*y + px[0]*py[3]*x - px[0]*py[3]*y + px[1]*py[0]*y - px[1]*py[0] + px[1]*py[2]*x - px[1]*py[3]*x + px[2]*py[0]*x - px[2]*py[0]*y - px[2]*py[1]*x - px[3]*py[0]*x + px[3]*py[0]*y + px[3]*py[1]*x);
Hx[0] += x1*x5;
Hx[1] += x4*x5;
Hx[2] += x2*x5;
Hx[3] += x3*x5;
			}

			UTOPIA_FUNCTION void gradient(
					// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Coefficients
				const T *UTOPIA_RESTRICT u,
				// Quadrature rule
				const T x,
				const T y,
				const T weight,
				T *UTOPIA_RESTRICT g) const
			{
				using namespace utopia::device;
			    // Automatically generated
				//FLOATING POINT OPS!
//	- Result: 4*ADDAUGMENTEDASSIGNMENT + 4*MUL
//	- Subexpressions: 10*ADD + 39*MUL + 11*SUB
T x0 = 1 - y;
T x1 = x0*(1 - x);
T x2 = x*y;
T x3 = 1 - x2;
T x4 = x*x0;
T x5 = weight*(u[0]*x1 + u[1]*x4 + u[2]*x2 + u[3]*x3)*(-px[0]*py[1]*y + px[0]*py[1] - px[0]*py[2]*x + px[0]*py[2]*y + px[0]*py[3]*x - px[0]*py[3]*y + px[1]*py[0]*y - px[1]*py[0] + px[1]*py[2]*x - px[1]*py[3]*x + px[2]*py[0]*x - px[2]*py[0]*y - px[2]*py[1]*x - px[3]*py[0]*x + px[3]*py[0]*y + px[3]*py[1]*x);
g[0] += x1*x5;
g[1] += x4*x5;
g[2] += x2*x5;
g[3] += x3*x5;
			}

			UTOPIA_FUNCTION void value(
					// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Coefficients
				const T *UTOPIA_RESTRICT u,
				// Quadrature rule
				const T x,
				const T y,
				const T weight,
				T &e
				) const
			{
				using namespace utopia::device;
			    // Automatically generated
				//TODO
			}

			UTOPIA_FUNCTION void eval(
					// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Coefficients
				const T *UTOPIA_RESTRICT u,
				// Quadrature rule
				const T x,
				const T y,
				const T weight,
				T &e,
				T *UTOPIA_RESTRICT g,
				T *UTOPIA_RESTRICT H) const
			{
				using namespace utopia::device;
			    // Automatically generated
				//FLOATING POINT OPS!
//	- Result: 20*ADDAUGMENTEDASSIGNMENT + 8*MUL + 3*POW
//	- Subexpressions: 10*ADD + 51*MUL + 2*POW + 11*SUB
T x0 = 1 - x;
T x1 = 1 - y;
T x2 = weight*(-px[0]*py[1]*y + px[0]*py[1] - px[0]*py[2]*x + px[0]*py[2]*y + px[0]*py[3]*x - px[0]*py[3]*y + px[1]*py[0]*y - px[1]*py[0] + px[1]*py[2]*x - px[1]*py[3]*x + px[2]*py[0]*x - px[2]*py[0]*y - px[2]*py[1]*x - px[3]*py[0]*x + px[3]*py[0]*y + px[3]*py[1]*x);
T x3 = pow(x1, 2)*x2;
T x4 = x*x0*x3;
T x5 = x*y;
T x6 = x0*x1;
T x7 = x2*x5*x6;
T x8 = 1 - x5;
T x9 = x2*x8;
T x10 = x6*x9;
T x11 = pow(x, 2);
T x12 = x11*x2;
T x13 = x1*x12*y;
T x14 = x*x1;
T x15 = x14*x9;
T x16 = x5*x9;
T x17 = u[0]*x6 + u[1]*x14 + u[2]*x5 + u[3]*x8;
T x18 = x17*x2;
H[0] += pow(x0, 2)*x3;
H[1] += x4;
H[2] += x7;
H[3] += x10;
H[4] += x4;
H[5] += x11*x3;
H[6] += x13;
H[7] += x15;
H[8] += x7;
H[9] += x13;
H[10] += x12*pow(y, 2);
H[11] += x16;
H[12] += x10;
H[13] += x15;
H[14] += x16;
H[15] += x2*pow(x8, 2);
g[0] += x18*x6;
g[1] += x14*x18;
g[2] += x18*x5;
g[3] += x17*x9;
			}

			//TODO

		};
	}

	namespace kokkos {
		template<class FE>
		using MassQuad4 = utopia::kokkos::AutoKernel<FE, utopia::kernels::Mass<utopia::kernels::Quad4<typename FE::Scalar, typename FE::Scalar>>, 2>;
	}
}

#endif // UTOPIA_TPL_MATERIAL_Mass_Quad4_2_IMPL_hpp
