#include <iostream>
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
//#include "utopia_FEDSLBaseSolverExamples.hpp"
#include "utopia_Biomechanics.hpp"
#include "utopia_UGMeshReader.hpp"

using namespace utopia;
using namespace std;
using namespace libMesh;


int main(const int argc, char *argv[]) 
{
	LibMeshInit init(argc, argv);
	

	if(argc > 3) {
		const std::string command = argv[1];

		if(command == "convert") {
			const std::string src = argv[2];
			const std::string dest = argv[3];

			Mesh mesh(init.comm());
			
			UGXMeshReader reader;
			if(!reader.read(src, mesh)) {
				return EXIT_FAILURE;
			}

			ExodusII_IO(mesh).write(dest);
			return EXIT_SUCCESS;
		}
	}

	// run_base_examples(init);
	// run_time_diff_examples(init);
	// run_mortar_examples(init);
	// run_least_squares_examples(init);
	// run_mixed_fe_space_example(init);
    //run_solver_ex(init);
//    run_biomechanics_example(init);
	return EXIT_SUCCESS;
}

