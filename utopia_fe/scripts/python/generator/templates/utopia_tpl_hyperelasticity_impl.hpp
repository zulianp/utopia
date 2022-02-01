#ifndef UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp

#include "utopia_hyperelasticity_{name}.hpp"

namespace utopia {{
	namespace kernels {{

		/** 
		 * Specialization of {name} for dimension {dim} 
		 */
		template<typename T>
		class {name}<T, {dim}> {{
		public:

			UTOPIA_FUNCTION static void hessian(const T mu,
				const T lmbda,
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *grad_trial,
				const T dx,
				T *UTOPIA_RESTRICT bf)
			{{
				using namespace utopia::device;
				// Automatically generated
				{hessian}
			}}

			UTOPIA_FUNCTION static void gradient(const T mu,
				const T lmbda,
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T dx,
				T *UTOPIA_RESTRICT lf)
			{{
				using namespace utopia::device;
			    // Automatically generated
				{gradient}
			}}

			UTOPIA_FUNCTION static void value(const T mu,
				const T lmbda,
				const T *UTOPIA_RESTRICT f,
				const T dx,
				T &e
				)
			{{
				using namespace utopia::device;
			    // Automatically generated
				{value}
			}}

			UTOPIA_FUNCTION void eval(const T mu,
				const T lmbda,
				const T *UTOPIA_RESTRICT f,
				const T *grad_test,
				const T *grad_trial,
				const T dx,
				T &e,
				T *UTOPIA_RESTRICT lf,
				T *UTOPIA_RESTRICT bf)
			{{
				using namespace utopia::device;
			    // Automatically generated
				{combined}
			}}

		}};;
	}}
}}

#endif // UTOPIA_TPL_HYPERELASTICITY_{name}_{dim}_IMPL_hpp
