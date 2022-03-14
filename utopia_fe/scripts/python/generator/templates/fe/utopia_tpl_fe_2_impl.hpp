#ifndef UTOPIA_TPL_FE_{name}_{dim}_IMPL_hpp
#define UTOPIA_TPL_FE_{name}_{dim}_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_{name}.hpp"

namespace utopia {{
	namespace kernels {{

		/**
		 * Specialization of {name} for dimension {dim}
		 */
		template<typename T, typename GeoT = T>
		class {name} {{
		public:
			static constexpr int Dim = {dim};
			using Result = typename utopia::MostDescriptive<T, GeoT>::Type;

			UTOPIA_FUNCTION static constexpr const char* class_name() {{ return "{name}"; }}

			UTOPIA_FUNCTION static constexpr Result measure(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				// Input quadrature point
				const T x,
				const T y)
			{{

				T measure_value;
				{measure}
				return measure_value;
			}}

			UTOPIA_FUNCTION static void gradient(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				// Input quadrature point
				const T x,
				const T y,
				// Output
				Result UTOPIA_RESTRICT*gx,
				Result UTOPIA_RESTRICT*gy)
			{{
				using namespace utopia::device;
				// Automatically generated
				{gradient}
			}}

			UTOPIA_FUNCTION static void value(
				const T x,
				const T y,
				Result UTOPIA_RESTRICT*f
				)
			{{
				using namespace utopia::device;
			    // Automatically generated
				{value}
			}}

		}};
	}}
}}

#endif // UTOPIA_TPL_FE_{name}_{dim}_IMPL_hpp
