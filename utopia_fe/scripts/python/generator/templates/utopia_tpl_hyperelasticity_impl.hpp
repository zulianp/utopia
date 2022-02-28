#ifndef UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_{name}.hpp"

namespace utopia {{
	namespace kernels {{

		/**
		 * Specialization of {name} for dimension {dim}
		 */
		template<typename T>
		class {name}<T, {dim}> {{
		public:
			static constexpr int Dim = {dim};

			UTOPIA_FUNCTION static constexpr const char* class_name() {{ return "{name}_{dim}"; }}

			class Params : public Configurable {{
			public:
				void read(Input &in) override
				{{
					{get_params}
				}}

				{fields}
			}};

			{name}(const Params &params)
			{{
				{set_params}
			}}

			UTOPIA_FUNCTION void hessian(
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *grad_trial,
				const T dx,
				T *UTOPIA_RESTRICT bf) const
			{{
				using namespace utopia::device;
				// Automatically generated
				{hessian}
			}}

			UTOPIA_FUNCTION void gradient(
				const T *UTOPIA_RESTRICT f,
				const T *UTOPIA_RESTRICT grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf) const
			{{
				using namespace utopia::device;
			    // Automatically generated
				{gradient}
			}}

			UTOPIA_FUNCTION void value(
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				) const
			{{
				using namespace utopia::device;
			    // Automatically generated
				{value}
			}}

			UTOPIA_FUNCTION void eval(
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *grad_trial,
				const T dx,
				T &e,
				T *UTOPIA_RESTRICT lf,
				T *UTOPIA_RESTRICT bf) const
			{{
				using namespace utopia::device;
			    // Automatically generated
				{combined}
			}}

			UTOPIA_FUNCTION void apply(
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *disp_grad,
				const T dx,
				T *UTOPIA_RESTRICT res) const
			{{
				using namespace utopia::device;
				// Automatically generated
				{apply_hessian}
			}}

			{fields}

		}};
	}}
}}

#endif // UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp
