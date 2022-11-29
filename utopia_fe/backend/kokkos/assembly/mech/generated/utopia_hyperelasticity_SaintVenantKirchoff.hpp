#ifndef UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_hpp
#define UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_kokkos_AutoHyperElasticity.hpp"


#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
	namespace kernels {

		template<typename T, int Dim>
		class SaintVenantKirchoff {
		public:
			static_assert(Dim < 0, "Automatically generated class -- SaintVenantKirchoff -- needs template specialization");

		};
	}


	namespace kokkos {
		template<class FunctionSpace, class FE, int Dim>
		using SaintVenantKirchoff = utopia::kokkos::AutoHyperElasticity<FunctionSpace, FE, utopia::kernels::SaintVenantKirchoff<typename FE::Scalar, Dim>>;
	}
}

#endif //UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_hpp
