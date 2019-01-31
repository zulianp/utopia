#include <iostream>

#include "utopia_libmesh.hpp"
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

#include "utopia_TimeDiffExamples.hpp"
#include "utopia_TimeDiffExamples.hpp"
#include "utopia_GeometryTest.hpp"
#include "utopia_UtopiaFETests.hpp"
#include "utopia_SemigeometricMultigridTest.hpp"
#include "utopia_SDCTest.hpp"
#include "utopia_MechTest.hpp"
#include "utopia_AssemblyTest.hpp"
#include "utopia_FSITest.hpp"
#include "utopia_MSHReaderTest.hpp"
#include "utopia_BoundaryIntegralTest.hpp"
#include "utopia_FormEvalTest.hpp"
#include "utopia_NonLinearElasticityTest.hpp"
#include "utopia_FEEvalTest.hpp"
#include "utopia_LeastSquaresHelmholtz.hpp"
#include "utopia_ContactApp.hpp"
#include "utopia_CoarsenerTest.hpp"
#include "utopia_EikonalEquationTest.hpp"
#include "utopia_TestVolume2SurfaceTransfer.hpp"
#include "utopia_VolumeInterpolationTest.hpp"
#include "utopia_WearSimulation.hpp"
#include "utopia_TransferApp.hpp"
#include "utopia_FractureFlowApp.hpp"
#include "utopia_RMTRApp.hpp"
#include "utopia_EnergyAssemblyTest.hpp"
#include "utopia_Intrepid2Test.hpp"
#include "utopia_IntersectTest.hpp"
#include "utopia_FETest.hpp"

#include "utopia_Grid2MeshTransferApp.hpp"

#include <functional>

#include "par_moonolith.hpp"

using namespace utopia;
using namespace std;
using namespace libMesh;


int main(int argc, char *argv[])
{
	//For debugginh with ddt
	MPI_Init(&argc, &argv);
	PETSC_COMM_WORLD = MPI_COMM_WORLD;

	Utopia::Init(argc, argv);
	MOONOLITH_PROFILING_BEGIN();

	{
		LibMeshInit init(argc, argv, PETSC_COMM_WORLD);

		std::map<std::string, std::function<void(LibMeshInit &)> > runners;
		// runners["base"] = run_base_examples;
		// runners["time_diff"] = run_time_diff_examples;
		// runners["least_squares"] = run_least_squares_examples;
		// runners["mixed_fe_space"] = run_mixed_fe_space_example;
	    // runners["biomechanics"] = run_biomechanics_example;
	    runners["geometry"] = run_geometry_test;
	    // runners["mortar"] = run_mortar_examples;
	    runners["tests"] = run_all_utopia_fe_tests;
	    runners["smg"] = run_semigeometric_multigrid_test;
	    runners["sdc"] = run_sdc_test;
	    // runners["wear"] = run_wear_test;
	    runners["mech"] = run_mech_test;
	    runners["asm"] = run_assembly_test;
	    runners["fsi"] = run_fsi_test;
	    runners["bit"] = run_boundary_integral_test;
	    runners["fet"] = run_form_eval_test;
	    runners["nle_test"] = run_non_linear_elasticity_test;
	    runners["test_msh_reader"] = test_msh_reader;
	    runners["fe_test"] = run_fe_eval_test;
	    runners["helm"] = run_form_least_squares_helmholtz;

	    // runners["ct"] = run_contact_test;
	    runners["coarsener_test"] = run_coarsener_test;
	    runners["eikonal"] = run_eikonal_equation_test;
	    runners["vol2surf"] = run_volume_to_surface_transfer_test;
	    runners["interp"] = run_volume_interpolation_test;
	    runners["energy"] = run_energy_test;
	    runners["intrepid2"] = run_intrepid2_test;
	    runners["isect"] = run_intersect_test;
	    runners["tests"] = run_all_tests;

	    // runners["coupled"] = run_coupled_equation_test;
	    //benchmarks
	    // runners["vt_benchmark"] = run_volume_transfer_benchmark;
	    // runners["vt_weak_scaling"] = run_weak_scaling_benchmark;


		for(int i = 1; i < argc; ++i) {
			const int ip1 = i+1;

			if(argv[i] == std::string("-verbose")) {
				Utopia::instance().set("verbose", "true");
				continue;
			}

			if(argv[i] == std::string("-r")) {
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
			} else if(argv[i] == std::string("-output_path")) {
				utopia::Utopia::instance().set("output_path", argv[ip1]);
				std::cout << "setting output_path to: " << argv[ip1] << std::endl;
			} else if(argv[i] == std::string("-data_path")) {
				utopia::Utopia::instance().set("data_path", argv[ip1]);
				std::cout << "setting data_path to: " << argv[ip1] << std::endl;
			} else if(argv[i] == std::string("-wear_sim")) {
				std::cout << argv[i] << " " << argv[ip1] << std::endl;
				//passing wear xml path to
				WearSimulation ws;
				ws.run(init, argv[ip1]);
			} else if(argv[i] == TransferApp::command()) {
				std::cout << argv[i] << " " << argv[ip1] << std::endl;

				TransferApp app;
				app.init(init);
				app.run(argv[ip1]);
			} else if(argv[i] == FractureFlowApp::command()) {
				std::cout << argv[i] << " " << argv[ip1] << std::endl;

				FractureFlowApp app;
				app.init(init);
				app.run(argv[ip1]);
			} else if(argv[i] == RMTRApp::command()) {
				std::cout << argv[i] << " " << argv[ip1] << std::endl;

				RMTRApp app;
				app.init(init);
				app.run(argv[ip1]);
			} else if(argv[i] == ContactApp::command()) {
				std::cout << argv[i] << " " << argv[ip1] << std::endl;

				ContactApp app;
				app.init(init);
				app.run(argv[ip1]);
			} else if(argv[i] == Grid2MeshTransferApp::command()) {
				std::cout << argv[i] << " " << argv[ip1] << std::endl;

				Grid2MeshTransferApp app;
				app.init(init);
				app.run(argv[ip1]);
			}

		}
	}

	MOONOLITH_PROFILING_END();
    return Utopia::Finalize();
}


