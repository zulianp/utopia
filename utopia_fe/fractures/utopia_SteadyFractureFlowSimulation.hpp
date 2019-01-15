#ifndef UTOPIA_STATIC_FRACTURE_FLOW_SIMULATION_HPP
#define UTOPIA_STATIC_FRACTURE_FLOW_SIMULATION_HPP

#include "libmesh/parallel_mesh.h"
#include "utopia_FractureFlow.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

	class SteadyFractureFlowSimulation : public Configurable {
	public:
		SteadyFractureFlowSimulation(libMesh::Parallel::Communicator &comm);
		void read(Input &is) override;
		void assemble_systems();
		void init_coupling();
		void write_output();
		bool solve();
		bool run();

		std::unique_ptr<FractureFlow> matrix;
		std::unique_ptr<FractureFlow> fracture_newtork;
		std::unique_ptr<FractureFlow> lagrange_multiplier;

		USparseMatrix A_m, A_f;
		UVector rhs_m, rhs_f;
		UVector x_m, x_f, lagr;

		USparseMatrix D, B, D_t, B_t;

		std::string solve_strategy;
		bool use_mg;
		int mg_sweeps;
		int mg_levels;
		bool plot_matrix;
		bool write_operators_to_disk;

	private:
		bool solve_cg_dual();
		bool solve_monolithic();
		bool solve_separate();
		bool solve_staggered();
	};
}

#endif //UTOPIA_STATIC_FRACTURE_FLOW_SIMULATION_HPP
