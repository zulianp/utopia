#ifndef UTOPIA_TPL_HYPERELASTICITY_GuccioneCosta_hpp
#define UTOPIA_TPL_HYPERELASTICITY_GuccioneCosta_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_kokkos_AutoHyperElasticity.hpp"


#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		template<typename T, int Dim>
		class GuccioneCosta {
		public:
			static_assert(Dim < 0, "Automatically generated class -- GuccioneCosta -- needs template specialization");

		};
	}

	namespace kokkos {
		template<class FunctionSpace, class FE, int Dim>
		using GuccioneCosta = utopia::kokkos::AutoHyperElasticity<FunctionSpace, FE, utopia::kernels::GuccioneCosta<typename FE::Scalar, Dim>>;
	}
}

#endif //UTOPIA_TPL_HYPERELASTICITY_GuccioneCosta_hpp
