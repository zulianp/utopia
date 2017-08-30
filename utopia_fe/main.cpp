#include <iostream>

#include "par_moonolith.hpp"
#include "utopia.hpp"

//fe extension
#include "utopia_fe.hpp"
#include <libmesh/const_function.h>
// #include <libmesh/vtk_io.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>

#include "MortarAssembler.hpp"
#include "utopia_Socket.hpp"

#include "utopia_FEDSLBaseExamples.hpp"
#include "utopia_FEDSLMortarExamples.hpp"
#include "utopia_FEDSLLeastSquaresExamples.hpp"
#include "utopia_MixedFESpaceExample.hpp"
#include "utopia_TimeDiffExamples.hpp"
#include "utopia_TimeDiffExamples.hpp"
#include "utopia_GeometryTest.hpp"
#include "utopia_Biomechanics.hpp"
#include "utopia_UtopiaFETests.hpp"

//#include "utopia_FEDSLBaseSolverExamples.hpp"
using namespace utopia;
using namespace std;
using namespace libMesh;


int main(const int argc, char *argv[]) 
{

	Utopia::Init(argc, argv);
	
	{
		LibMeshInit init(argc, argv, PETSC_COMM_WORLD);
		run_all_utopia_fe_tests(init);
	}
	
	// run_base_examples(init);
	// run_time_diff_examples(init);
	
	// run_least_squares_examples(init);
	// run_mixed_fe_space_example(init);
    //run_solver_ex(init);
    // run_biomechanics_example(init);
    // run_geometry_test(init);
    // run_mortar_examples(init);
    return Utopia::Finalize();
	// return EXIT_SUCCESS;
}

