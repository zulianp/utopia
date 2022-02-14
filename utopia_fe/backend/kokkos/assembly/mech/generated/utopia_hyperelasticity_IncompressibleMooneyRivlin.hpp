#ifndef UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_hpp
#define UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_kokkos_MixedHyperElasticity.hpp"


#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		template<typename T, int Dim>
		class IncompressibleMooneyRivlin {
		public:
			static_assert(Dim < 0, "Automatically generated class -- IncompressibleMooneyRivlin -- needs template specialization");

		};
	}


	namespace kokkos {
		template<class FE, int Dim>
		using IncompressibleMooneyRivlin = utopia::kokkos::MixedHyperElasticity<FE, utopia::kernels::IncompressibleMooneyRivlin<typename FE::Scalar, Dim>>;
	}
}

#endif //UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_hpp
