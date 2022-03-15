#ifndef UTOPIA_TPL_FE_AxisAlignedQuad4_2_IMPL_hpp
#define UTOPIA_TPL_FE_AxisAlignedQuad4_2_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"

// #include "utopia_fe_AxisAlignedQuad4.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of AxisAlignedQuad4 for dimension 2
		 */
		template<typename T, typename GeoT = T>
		class AxisAlignedQuad4 {
		public:
			static constexpr int Dim = 2;
			static constexpr int NNodes = 4;
			static constexpr int Order = 1;

			using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

			UTOPIA_FUNCTION static constexpr const char* class_name() { return "AxisAlignedQuad4"; }

			UTOPIA_INLINE_FUNCTION static constexpr int dim() 
			{
				return Dim;
			}

			UTOPIA_INLINE_FUNCTION static constexpr int n_nodes() 
			{
				return NNodes;
			}

			UTOPIA_INLINE_FUNCTION static constexpr int order() 
			{
				return Order;
			}

			UTOPIA_FUNCTION static constexpr Result measure(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Input quadrature point
				const T x,
				const T y)
			{
				T measure_value;
				//FLOATING POINT OPS!
//	- Result: ADD + ASSIGNMENT + 4*MUL
//	- Subexpressions: 0
measure_value = px[0]*py[0] - px[0]*py[2] - px[2]*py[0] + px[2]*py[2];
				return measure_value;
			}

			UTOPIA_FUNCTION static void jacobian(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Input quadrature point
				const T x,
				const T y,
				GeoT *UTOPIA_RESTRICT J)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 2*ADD + 4*ASSIGNMENT + 2*MUL
//	- Subexpressions: 0
J[0] = -px[0] + px[2];
J[1] = 0;
J[2] = 0;
J[3] = -py[0] + py[2];
			}

			UTOPIA_FUNCTION static void jacobian_inverse(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Input quadrature point
				const T x,
				const T y,
				GeoT *UTOPIA_RESTRICT J_inv)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 2*ADD + 4*ASSIGNMENT + 4*MUL + 2*POW
//	- Subexpressions: 0
J_inv[0] = -1/(px[0] - px[2]);
J_inv[1] = 0;
J_inv[2] = 0;
J_inv[3] = -1/(py[0] - py[2]);
			}

			UTOPIA_FUNCTION static void transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Input quadrature point
				const T x,
				const T y,
				GeoT &tx,
				GeoT &ty)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 4*ADD + 2*ASSIGNMENT + 4*MUL
//	- Subexpressions: 0
tx = px[0] - x*(px[0] - px[2]);
ty = py[0] - y*(py[0] - py[2]);
			}

			UTOPIA_FUNCTION static void inverse_transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Input quadrature point
				const T tx,
				const T ty,
				GeoT &x,
				GeoT &y)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 4*ADD + 2*ASSIGNMENT + 6*MUL + 2*POW
//	- Subexpressions: 0
x = (px[0] - tx)/(px[0] - px[2]);
y = (py[0] - ty)/(py[0] - py[2]);
			}

			UTOPIA_FUNCTION static void gradient(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				// Input quadrature point
				const T x,
				const T y,
				// Output
				Result *UTOPIA_RESTRICT gx,
				Result *UTOPIA_RESTRICT gy)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: ADD + 8*ASSIGNMENT + 6*MUL
//	- Subexpressions: 2*DIV + 2*MUL + 3*SUB
T x0 = y - 1;
T x1 = 1.0/(px[0] - px[2]);
T x2 = 1.0/(py[0] - py[2]);
T x3 = x*x2;
T x4 = x1*y;
gx[0] = -x0*x1;
gy[0] = x2*(1 - x);
gx[1] = x0*x1;
gy[1] = x3;
gx[2] = -x4;
gy[2] = -x3;
gx[3] = x4;
gy[3] = x3;
			}

			UTOPIA_FUNCTION static void value(
				const T x,
				const T y,
				Result *UTOPIA_RESTRICT f
				)
			{
				using namespace utopia::device;
			    // Automatically generated
				//FLOATING POINT OPS!
//	- Result: 2*ADD + 4*ASSIGNMENT + 4*MUL
//	- Subexpressions: MUL + SUB
T x0 = 1 - y;
T x1 = x*y;
f[0] = x0*(1 - x);
f[1] = x*x0;
f[2] = x1;
f[3] = 1 - x1;
			}


		UTOPIA_FUNCTION static void eval(
			// Element coordinates
			const GeoT *UTOPIA_RESTRICT px,
			const GeoT *UTOPIA_RESTRICT py,
			// Input quadrature point
			const T x,
			const T y,
			// Output
			Result *UTOPIA_RESTRICT f,
			Result *UTOPIA_RESTRICT gx,
			Result *UTOPIA_RESTRICT gy,
			T &measure_value)
		{
			using namespace utopia::device;
			// Automatically generated
			//FLOATING POINT OPS!
//	- Result: 2*ADD + 13*ASSIGNMENT + 12*MUL
//	- Subexpressions: 2*DIV + 3*MUL + NEG + 4*SUB
T x0 = 1 - x;
T x1 = y - 1;
T x2 = -x1;
T x3 = x*y;
T x4 = 1.0/(px[0] - px[2]);
T x5 = 1.0/(py[0] - py[2]);
T x6 = x*x5;
T x7 = x4*y;
f[0] = x0*x2;
f[1] = x*x2;
f[2] = x3;
f[3] = 1 - x3;
measure_value = px[0]*py[0] - px[0]*py[2] - px[2]*py[0] + px[2]*py[2];
gx[0] = x2*x4;
gy[0] = x0*x5;
gx[1] = x1*x4;
gy[1] = x6;
gx[2] = -x7;
gy[2] = -x6;
gx[3] = x7;
gy[3] = x6;
		}

		};
	}
}

#endif // UTOPIA_TPL_FE_AxisAlignedQuad4_2_IMPL_hpp
