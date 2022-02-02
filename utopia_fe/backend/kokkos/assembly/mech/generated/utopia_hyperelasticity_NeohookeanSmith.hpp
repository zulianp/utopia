#ifndef UTOPIA_TPL_HYPERELASTICITY_NeohookeanSmith_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeohookeanSmith_hpp

#include "utopia_Algorithms.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		template<typename T, int Dim>
		class NeohookeanSmith {
		public:
			static_assert(Dim < 0, "Automatically generated class -- NeohookeanSmith -- needs template specialization");

		};
	}
}

#endif //UTOPIA_TPL_HYPERELASTICITY_NeohookeanSmith_hpp
