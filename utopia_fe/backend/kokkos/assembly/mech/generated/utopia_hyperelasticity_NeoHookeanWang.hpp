#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_kokkos_AutoHyperElasticity.hpp"


#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		template<typename T, int Dim>
		class NeoHookeanWang {
		public:
			static_assert(Dim < 0, "Automatically generated class -- NeoHookeanWang -- needs template specialization");

		};
	}


	namespace kokkos {
		template<class FunctionSpace, class FE, int Dim>
		using NeoHookeanWang = utopia::kokkos::AutoHyperElasticity<FunctionSpace, FE, utopia::kernels::NeoHookeanWang<typename FE::Scalar, Dim>>;
	}
}

#endif //UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_hpp
