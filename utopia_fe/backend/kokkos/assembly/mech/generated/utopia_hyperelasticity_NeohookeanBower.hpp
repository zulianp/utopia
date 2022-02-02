#ifndef UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_hpp

#include "utopia_Algorithms.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		template<typename T, int Dim>
		class NeohookeanBower {
		public:
			static_assert(Dim < 0, "Automatically generated class -- NeohookeanBower -- needs template specialization");

		};
	}
}

#endif //UTOPIA_TPL_HYPERELASTICITY_NeohookeanBower_hpp
