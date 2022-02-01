#ifndef UTOPIA_TPL_HYPERELASTICITY_{name}_hpp
#define UTOPIA_TPL_HYPERELASTICITY_{name}_hpp

#include "utopia_Algorithms.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {{
	namespace kernels {{

		template<typename T, int Dim>
		class {name} {{
		public:
			static_assert(Dim < 0, "Automatically generated class -- {name} -- needs template specialization");

		}};
	}}
}}

#endif //UTOPIA_TPL_HYPERELASTICITY_{name}_hpp
