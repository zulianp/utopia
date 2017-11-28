#include "utopia_UtopiaFETests.hpp"
#include "utopia_SemigeometricMultigridTest.hpp"
#include "utopia_ContactTest.hpp"
#include "utopia_FormEvalTest.hpp"
#include "utopia_LibMeshBackendTest.hpp"
#include "utopia_FSITest.hpp"

namespace utopia {
	void run_all_utopia_fe_tests(libMesh::LibMeshInit &init)
	{
		// run_libmesh_backend_test(init);
		// run_contact_test(init);
		run_form_eval_test(init);
		// run_fsi_test(init);
		
		// run_semigeometric_multigrid_test(init);
	}
}

