#include "utopia_UtopiaFETests.hpp"
#include "utopia_SemigeometricMultigridTest.hpp"
#include "utopia_ContactTest.hpp"

namespace utopia {
	void run_all_utopia_fe_tests(libMesh::LibMeshInit &init)
	{
		//run_contact_test(init);
		 run_semigeometric_multigrid_test(init);
	}
}
