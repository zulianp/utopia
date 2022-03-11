#ifndef UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp

#include "utopia_Input.hpp"

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

			UTOPIA_FUNCTION static constexpr const char* class_name() {{ return "{name}"; }}

			UTOPIA_FUNCTION static constexpr T measure(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				// Input quadrature point
				const T x,
				const T y)
			{
				{{measure}}
			}

			UTOPIA_FUNCTION static void gradient(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				// Input quadrature point
				const T x,
				const T y,
				// Output
				T UTOPIA_RESTRICT*gx,
				T UTOPIA_RESTRICT*gy)
			{{
				using namespace utopia::device;
				// Automatically generated
				{hessian}
			}}

			UTOPIA_FUNCTION static void value(
				const T x,
				const T y,
				T UTOPIA_RESTRICT*f,
				)
			{{
				using namespace utopia::device;
			    // Automatically generated
				{value}
			}}

		}};
	}}
}}

#endif // UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp
