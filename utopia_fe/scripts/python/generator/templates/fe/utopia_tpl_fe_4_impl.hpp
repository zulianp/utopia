#ifndef UTOPIA_TPL_FE_{name}_{dim}_IMPL_hpp
#define UTOPIA_TPL_FE_{name}_{dim}_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Algorithms.hpp"

// #include "utopia_fe_{name}.hpp"

#include <cassert>

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {{
	namespace kernels {{

		/**
		 * Specialization of {name} for dimension {dim}
		 */
		template<typename T, typename GeoT>
		class {name} {{
		public:
			static constexpr int Dim = {dim};
			static constexpr int NNodes = {nnodes};
			static constexpr int Order = {order};

			using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

			UTOPIA_FUNCTION static constexpr const char* class_name() {{ return "{name}"; }}

			UTOPIA_INLINE_FUNCTION static constexpr int dim() 
			{{
				return Dim;
			}}

			UTOPIA_INLINE_FUNCTION static constexpr int n_nodes() 
			{{
				return NNodes;
			}}

			UTOPIA_INLINE_FUNCTION static constexpr int order() 
			{{
				return Order;
			}}

			UTOPIA_FUNCTION static Result measure(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t)
			{{
				T measure_value = 0;
				{measure}
				return measure_value;
			}}

			UTOPIA_FUNCTION static void jacobian(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				GeoT *UTOPIA_RESTRICT J)
			{{
				using namespace utopia::device;
				// Automatically generated
				{jacobian}
			}}

			UTOPIA_FUNCTION static void jacobian_inverse(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				GeoT *UTOPIA_RESTRICT J_inv)
			{{
				using namespace utopia::device;
				// Automatically generated
				{jacobian_inverse}
			}}

			UTOPIA_FUNCTION static void transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				GeoT &tx,
				GeoT &ty,
				GeoT &tz,
				GeoT &tt)
			{{
				using namespace utopia::device;
				// Automatically generated
				{transform}
			}}

			UTOPIA_FUNCTION static void inverse_transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T tx,
				const T ty,
				const T tz,
				const T tt,
				GeoT &x,
				GeoT &y,
				GeoT &z,
				GeoT &t)
			{{
				using namespace utopia::device;
				// Automatically generated
				{inverse_transform}
			}}

			UTOPIA_FUNCTION static void gradient(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				// Output
				Result *UTOPIA_RESTRICT gx,
				Result *UTOPIA_RESTRICT gy,
				Result *UTOPIA_RESTRICT gz,
				Result *UTOPIA_RESTRICT gt)
			{{
				using namespace utopia::device;
				// Automatically generated
				{gradient}
			}}

			UTOPIA_FUNCTION static void value(
				const T x,
				const T y,
				const T z,
				const T t,
				Result *UTOPIA_RESTRICT f
				)
			{{
				using namespace utopia::device;
			    // Automatically generated
				{value}
			}}

			UTOPIA_FUNCTION static void eval(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				const GeoT *UTOPIA_RESTRICT py,
				const GeoT *UTOPIA_RESTRICT pz,
				const GeoT *UTOPIA_RESTRICT pt,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const T t,
				// Output
				Result *UTOPIA_RESTRICT f,
				Result *UTOPIA_RESTRICT gx,
				Result *UTOPIA_RESTRICT gy,
				Result *UTOPIA_RESTRICT gz,
				Result *UTOPIA_RESTRICT gt,
				T &measure_value)
			{{
				using namespace utopia::device;
				// Automatically generated
				{combined}
			}}

		}};
	}}
}}

#endif // UTOPIA_TPL_FE_{name}_{dim}_IMPL_hpp
