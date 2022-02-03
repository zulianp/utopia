#ifndef UTOPIA_TPL_HYPERELASTICITY_Fung_hpp
#define UTOPIA_TPL_HYPERELASTICITY_Fung_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_kokkos_AutoHyperElasticity.hpp"


#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		template<typename T, int Dim>
		class Fung {
		public:
			static_assert(Dim < 0, "Automatically generated class -- Fung -- needs template specialization");

		};
	}


	namespace kokkos {
		template<class FE, int Dim>
		using Fung = utopia::kokkos::AutoHyperElasticity<FE, utopia::kernels::Fung<typename FE::Scalar, Dim>>;
	}
}

#endif //UTOPIA_TPL_HYPERELASTICITY_Fung_hpp
