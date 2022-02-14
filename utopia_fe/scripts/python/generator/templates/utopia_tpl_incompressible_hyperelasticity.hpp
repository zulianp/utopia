#ifndef UTOPIA_TPL_HYPERELASTICITY_{name}_hpp
#define UTOPIA_TPL_HYPERELASTICITY_{name}_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_kokkos_MixedHyperElasticity.hpp"


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


	namespace kokkos {{
		template<class FE, int Dim>
		using {name} = utopia::kokkos::MixedHyperElasticity<FE, utopia::kernels::{name}<typename FE::Scalar, Dim>>;
	}}
}}

#endif //UTOPIA_TPL_HYPERELASTICITY_{name}_hpp
