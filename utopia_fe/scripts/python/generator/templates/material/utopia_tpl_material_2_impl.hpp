#ifndef UTOPIA_TPL_MATERIAL_{name}_{dim}_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_{name}_{dim}_IMPL_hpp

#include "utopia_Input.hpp"
#include "utopia_Algorithms.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_{trial}_{dim}.hpp"
#include "utopia_material_{name}.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif


namespace utopia {{
	namespace kernels {{

		/**
		 * Specialization of {name} for symmetric element pair trial=test={trial}
		 */
		template<typename T, typename GeoT>
		class {name}<{trial}<T, GeoT>> {{
		public:
			using ElemT = {trial}<T, GeoT>;
			static constexpr int Dim = ElemT::Dim;

			UTOPIA_FUNCTION static constexpr const char* class_name() {{ return "{name}<{trial}>"; }}

			class Params : public Configurable {{
			public:
				void read(Input &in) override
				{{
					{get_params}
				}}

				{fields}
			}};

			{name}(const Params &params = Params())
			{{
				{set_params}
			}}

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
			{{
				using namespace utopia::device;
				// Automatically generated
				{hessian}
			}}

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
			{{
				using namespace utopia::device;
				// Automatically generated
				{apply_hessian}
			}}

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
			{{
				using namespace utopia::device;
			    // Automatically generated
				{gradient}
			}}

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
			{{
				using namespace utopia::device;
			    // Automatically generated
				{value}
			}}

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
			{{
				using namespace utopia::device;
			    // Automatically generated
				{combined}
			}}

			{fields}

		}};
	}}

	namespace kokkos {{
		template<class FE>
		using {name}{trial} = utopia::kokkos::AutoKernel<FE, utopia::kernels::{name}<utopia::kernels::{trial}<typename FE::Scalar, typename FE::Scalar>>, {dim}>;
	}}
}}

#endif // UTOPIA_TPL_MATERIAL_{name}_{dim}_IMPL_hpp
