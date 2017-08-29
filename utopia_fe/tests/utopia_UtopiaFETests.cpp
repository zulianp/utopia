#include "utopia_UtopiaFETests.hpp"
#include "utopia_SemigeometricMultigridTest.hpp"

namespace utopia {
	void run_all_utopia_fe_tests(libMesh::LibMeshInit &init)
	{
		run_semigeometric_multigrid_test(init);
	}
}
