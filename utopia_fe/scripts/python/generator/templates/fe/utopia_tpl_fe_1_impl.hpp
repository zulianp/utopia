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

			UTOPIA_INLINE_FUNCTION static constexpr T reference_measure()
			{{
				return {reference_measure};
			}}

			UTOPIA_FUNCTION static constexpr Result measure(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x)
			{{
				T measure_value;
				{measure}
				return measure_value;
			}}

			UTOPIA_FUNCTION static void jacobian(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x,
				GeoT *UTOPIA_RESTRICT J)
			{{
				using namespace utopia::device;
				// Automatically generated
				{jacobian}
			}}

			UTOPIA_FUNCTION static void jacobian_inverse(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x,
				GeoT *UTOPIA_RESTRICT J_inv)
			{{
				using namespace utopia::device;
				// Automatically generated
				{jacobian_inverse}
			}}

			UTOPIA_FUNCTION static void transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x,
				GeoT &tx)
			{{
				using namespace utopia::device;
				// Automatically generated
				{transform}
			}}

			UTOPIA_FUNCTION static void inverse_transform(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T tx,
				GeoT &x)
			{{
				using namespace utopia::device;
				// Automatically generated
				{inverse_transform}
			}}

			UTOPIA_FUNCTION static void gradient(
				// Element coordinates
				const GeoT *UTOPIA_RESTRICT px,
				// Input quadrature point
				const T x,
				// Output
				Result *UTOPIA_RESTRICT gx)
			{{
				using namespace utopia::device;
				// Automatically generated
				{gradient}
			}}

			UTOPIA_FUNCTION static void value(
				const T x,
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
			// Input quadrature point
			const T x,
			// Output
			Result *UTOPIA_RESTRICT f,
			Result *UTOPIA_RESTRICT gx,
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
