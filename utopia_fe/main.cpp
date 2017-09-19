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
#include <functional>

//#include "utopia_FEDSLBaseSolverExamples.hpp"
using namespace utopia;
using namespace std;
using namespace libMesh;


int main(const int argc, char *argv[]) 
{

	Utopia::Init(argc, argv);
	
	{
		LibMeshInit init(argc, argv, PETSC_COMM_WORLD);

		std::map<std::string, std::function<void(LibMeshInit &)> > runners;
		runners["base"] = run_base_examples;
		runners["time_diff"] = run_time_diff_examples;
		runners["least_squares"] = run_least_squares_examples;
		runners["mixed_fe_space"] = run_mixed_fe_space_example;
	    runners["biomechanics"] = run_biomechanics_example;
	    runners["geometry"] = run_geometry_test;
	    runners["mortar"] = run_mortar_examples;
	    runners["tests"] = run_all_utopia_fe_tests;

		for(int i = 1; i < argc; ++i) {
			if(argv[i] == std::string("-r")) {
				const int ip1 = i+1;
				if(ip1 < argc) {
					auto it = runners.find(argv[ip1]);
					if(it == runners.end()) {
						std::cerr << "[Error] " << argv[ip1] << " not found" << std::endl;
					} else {
						std::cout << "--------------------------------------------" << std::endl;
						std::cout << "[Status] Running: " << argv[ip1] << std::endl;
						std::cout << "--------------------------------------------" << std::endl;
						it->second(init);
						std::cout << "--------------------------------------------" << std::endl;
						std::cout << "--------------------------------------------" << std::endl;
					}
				} else {
					std::cerr << "[Error] run requires an input string" << std::endl;
				}
			} else if(argv[i] == std::string("-h")) {
				std::cout << "--------------------------------------------" << std::endl;
				std::cout << "--------------------------------------------" << std::endl;
				std::cout << "-r <runner name>\n";
				std::cout << "Available runners:\n";

				for(auto &r : runners) {
					std::cout << "\t" << r.first << std::endl;
				}

				std::cout << "--------------------------------------------" << std::endl;
				std::cout << "--------------------------------------------" << std::endl;
			}
		}
	}

    return Utopia::Finalize();
}

