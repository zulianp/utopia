#ifndef UTOPIA_TPL_FE_Line2_1_IMPL_hpp
#define UTOPIA_TPL_FE_Line2_1_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"

// #include "utopia_fe_Line2.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		/**
		 * Specialization of Line2 for dimension 1
		 */
		template<typename T, typename GeoT>
		class Line2 {
		public:
			static constexpr int Dim = 1;
			static constexpr int NNodes = 2;
			static constexpr int Order = 1;

			using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

			UTOPIA_FUNCTION static constexpr const char* class_name() { return "Line2"; }

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

			UTOPIA_INLINE_FUNCTION static constexpr T reference_measure()
			{
				return 2;
			}

			UTOPIA_FUNCTION static constexpr Result measure(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x)
			{
				T measure_value;
				
// Unused variables
UTOPIA_UNUSED(x);
//FLOATING POINT OPS!
//	- Result: ADD + ASSIGNMENT + 2*MUL
//	- Subexpressions: 0
measure_value = -1.0/2.0*px[0] + (1.0/2.0)*px[1];
				return measure_value;
			}

			UTOPIA_FUNCTION static void jacobian(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x,
				GeoT *UTOPIA_RESTRICT J)
			{
				using namespace utopia::device;
				// Automatically generated
				
// Unused variables
UTOPIA_UNUSED(x);
//FLOATING POINT OPS!
//	- Result: ADD + ASSIGNMENT + 2*MUL
//	- Subexpressions: 0
J[0] = -1.0/2.0*px[0] + (1.0/2.0)*px[1];
			}

			UTOPIA_FUNCTION static void jacobian_inverse(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x,
				GeoT *UTOPIA_RESTRICT J_inv)
			{
				using namespace utopia::device;
				// Automatically generated
				
// Unused variables
UTOPIA_UNUSED(x);
//FLOATING POINT OPS!
//	- Result: ADD + ASSIGNMENT + 2*MUL + POW
//	- Subexpressions: 0
J_inv[0] = 1.0/(-1.0/2.0*px[0] + (1.0/2.0)*px[1]);
			}

			UTOPIA_FUNCTION static void transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x,
				GeoT &tx)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 3*ADD + ASSIGNMENT + 3*MUL
//	- Subexpressions: DIV
T x0 = (1.0/2.0)*x;
tx = px[0]*(0.5 - x0) + px[1]*(x0 + 0.5);
			}

			UTOPIA_FUNCTION static void inverse_transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T tx,
				GeoT &x)
			{
				using namespace utopia::device;
				// Automatically generated
				//FLOATING POINT OPS!
//	- Result: 2*ADD + ASSIGNMENT + 3*MUL + POW
//	- Subexpressions: 0
x = (px[0] + px[1] - 2.0*tx)/(px[0] - px[1]);
			}

			UTOPIA_FUNCTION static void gradient(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x,
				// Output
				Result *UTOPIA_RESTRICT gx)
			{
				using namespace utopia::device;
				// Automatically generated
				
// Unused variables
UTOPIA_UNUSED(x);
//FLOATING POINT OPS!
//	- Result: 2*ASSIGNMENT + MUL
//	- Subexpressions: DIV + SUB
T x0 = (1.0/2.0)/(-1.0/2.0*px[0] + (1.0/2.0)*px[1]);
gx[0] = -x0;
gx[1] = x0;
			}

			UTOPIA_FUNCTION static void value(
				const T x,
				Result *UTOPIA_RESTRICT f
				)
			{
				using namespace utopia::device;
			    // Automatically generated
				//FLOATING POINT OPS!
//	- Result: 2*ADD + 2*ASSIGNMENT + MUL
//	- Subexpressions: DIV
T x0 = (1.0/2.0)*x;
f[0] = 0.5 - x0;
f[1] = x0 + 0.5;
			}


		UTOPIA_FUNCTION static void eval(
			// Element coordinates
			const GeoT *UTOPIA_RESTRICT px,
			// Input quadrature point
			const T x,
			// Output
			Result *UTOPIA_RESTRICT f,
			Result *UTOPIA_RESTRICT gx,
			T &measure_value)
		{
			using namespace utopia::device;
			// Automatically generated
			//FLOATING POINT OPS!
//	- Result: 2*ADD + 5*ASSIGNMENT + 2*MUL
//	- Subexpressions: 4*DIV + MUL + SUB
T x0 = (1.0/2.0)*x;
T x1 = -1.0/2.0*px[0] + (1.0/2.0)*px[1];
T x2 = (1.0/2.0)/x1;
f[0] = 0.5 - x0;
f[1] = x0 + 0.5;
measure_value = x1;
gx[0] = -x2;
gx[1] = x2;
		}

		};
	}
}

#endif // UTOPIA_TPL_FE_Line2_1_IMPL_hpp
