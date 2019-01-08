#ifndef UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP
#define UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP

#include "utopia_SteadyFractureFlowSimulation.hpp"

namespace utopia {

	class FractureFlowTransportSimulation : public Configurable {
	public:
		FractureFlowTransportSimulation(libMesh::Parallel::Communicator &comm);
		void read(utopia::Input &in) override;
		void run();
		void compute_velocity();
		void append_aux_systems();
		void write_output();

	private:
		class Velocity {
		public:
			void init(const UVector &pressure, FractureFlow &flow);
			void update_output();

			std::unique_ptr<ProductFunctionSpace<LibMeshFunctionSpace>> space;
			UVector velocity;
			bool lump_mass_matrix;
		};

		SteadyFractureFlowSimulation steady_flow_;
		USparseMatrix upwind_f_;

		Velocity velocity_m_;
		Velocity velocity_f_;


	};

}

#endif //UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP
