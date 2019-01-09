#ifndef UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP
#define UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP

#include "utopia_SteadyFractureFlowSimulation.hpp"

namespace utopia {

	class FractureFlowTransportSimulation : public Configurable {
	public:
		FractureFlowTransportSimulation(libMesh::Parallel::Communicator &comm);
		void read(utopia::Input &in) override;
		void run();
		void compute_transport();
		void append_aux_systems();
		void write_output();

	private:
		class Transport {
		public:
			void init(const UVector &pressure, FractureFlow &flow);
			void update_output();
			void assemble_system();

			std::unique_ptr<ProductFunctionSpace<LibMeshFunctionSpace>> space;
			UVector velocity;
			bool lump_mass_matrix;
			USparseMatrix mass_matrix;
			USparseMatrix gradient_matrix;
			USparseMatrix system_matrix;
			UVector mass_vector;
			double dt;
		};

		SteadyFractureFlowSimulation steady_flow_;
		USparseMatrix upwind_f_;

		Transport transport_m_;
		Transport transport_f_;
		double simulation_time_;

	};

}

#endif //UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP
