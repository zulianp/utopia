#ifndef UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_hpp
#define UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_kokkos_AutoHyperElasticity.hpp"


#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		template<typename T, int Dim>
		class MooneyRivlin {
		public:
			static_assert(Dim < 0, "Automatically generated class -- MooneyRivlin -- needs template specialization");

		};
	}


	namespace kokkos {
		template<class FE, int Dim>
		using MooneyRivlin = utopia::kokkos::AutoHyperElasticity<FE, utopia::kernels::MooneyRivlin<typename FE::Scalar, Dim>>;
	}
}

#endif //UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_hpp
