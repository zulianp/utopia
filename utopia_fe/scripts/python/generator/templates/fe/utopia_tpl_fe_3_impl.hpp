#ifndef UTOPIA_TPL_FE_{name}_{dim}_IMPL_hpp
#define UTOPIA_TPL_FE_{name}_{dim}_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_{name}.hpp"

#include <cassert>

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
				const GeoT UTOPIA_RESTRICT*pz,
				// Input quadrature point
				const T x,
				const T y,
				const T z)
			{{
				T measure_value;
				{measure}
				return measure_value;
			}}

			UTOPIA_FUNCTION static void jacobian(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				const GeoT UTOPIA_RESTRICT*pz,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const GeoT UTOPIA_RESTRICT*J)
			{{
				using namespace utopia::device;
				// Automatically generated
				{jacobian}
			}}

			UTOPIA_FUNCTION static void jacobian_inverse(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				const GeoT UTOPIA_RESTRICT*pz,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const GeoT UTOPIA_RESTRICT*J_inv)
			{{
				using namespace utopia::device;
				// Automatically generated
				{jacobian_inverse}
			}}

			UTOPIA_FUNCTION static void transform(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				const GeoT UTOPIA_RESTRICT*pz,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const GeoT UTOPIA_RESTRICT*tx,
				const GeoT UTOPIA_RESTRICT*ty,
				const GeoT UTOPIA_RESTRICT*tz)
			{{
				using namespace utopia::device;
				// Automatically generated
				{transform}
			}}

			UTOPIA_FUNCTION static void inverse_transform(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				const GeoT UTOPIA_RESTRICT*pz,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				const GeoT UTOPIA_RESTRICT*tx,
				const GeoT UTOPIA_RESTRICT*ty,
				const GeoT UTOPIA_RESTRICT*tz)
			{{
				using namespace utopia::device;
				// Automatically generated
				{inverse_transform}
			}}

			UTOPIA_FUNCTION static void gradient(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				const GeoT UTOPIA_RESTRICT*pz,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				// Output
				Result UTOPIA_RESTRICT*gx,
				Result UTOPIA_RESTRICT*gy,
				Result UTOPIA_RESTRICT*gz)
			{{
				using namespace utopia::device;
				// Automatically generated
				{gradient}
			}}

			UTOPIA_FUNCTION static void value(
				const T x,
				const T y,
				const T z,
				Result UTOPIA_RESTRICT*f
				)
			{{
				using namespace utopia::device;
			    // Automatically generated
				{value}
			}}


			UTOPIA_FUNCTION static void eval(
				// Element coordinates
				const GeoT UTOPIA_RESTRICT*px,
				const GeoT UTOPIA_RESTRICT*py,
				const GeoT UTOPIA_RESTRICT*pz,
				// Input quadrature point
				const T x,
				const T y,
				const T z,
				// Output
				Result UTOPIA_RESTRICT*f,
				Result UTOPIA_RESTRICT*gx,
				Result UTOPIA_RESTRICT*gy,
				Result UTOPIA_RESTRICT*gz,
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
