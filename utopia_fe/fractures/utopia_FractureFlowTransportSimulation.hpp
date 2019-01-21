#ifndef UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP
#define UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP

#include "utopia_SteadyFractureFlowSimulation.hpp"
#include "utopia_UIForcingFunction.hpp"

namespace utopia {

	class FractureFlowTransportSimulation : public Configurable {
	public:
		FractureFlowTransportSimulation(libMesh::Parallel::Communicator &comm);
		void read(utopia::Input &in) override;
		void run();
		void compute_transport();
		void compute_transport_separate();
		void compute_transport_monolithic();
		void append_aux_systems();
		void write_output();
		void compute_upwind_operator();

	private:
		class Transport : public Configurable {
		public:
			inline void set_steady_state_function_space(UIFunctionSpace<LibMeshFunctionSpace> &V)
			{
				steady_state_function_space = make_ref(V);
				space = utopia::make_unique<UIFunctionSpace<LibMeshFunctionSpace>>(V.mesh(), V.subspace(0).equation_systems_ptr());
			}

			Transport();

			void init(const UVector &pressure, FractureFlow &flow);
			void update_output();
			void assemble_system(FractureFlow &flow);
			void remove_mass(const UVector &in, UVector &out);
			void add_mass(const UVector &in, UVector &out) const;
			void read(Input &in) override;
			void constrain_concentration(UVector &vec);
			void assemble_aux_quantities(FractureFlow &flow);
			
			void post_process_time_step(FractureFlow &flow);

			std::shared_ptr<UIFunctionSpace<LibMeshFunctionSpace>> steady_state_function_space;
			std::unique_ptr<UIFunctionSpace<LibMeshFunctionSpace>> space;
			std::unique_ptr<UIForcingFunction<LibMeshFunctionSpace, UVector>> forcing_function;
			ProductFunctionSpace<LibMeshFunctionSpace> aux_space;
			
			UVector concentration;
			bool lump_mass_matrix;
			bool h1_regularization;
			bool use_upwinding;
			double regularization_parameter;
			USparseMatrix mass_matrix;
			USparseMatrix gradient_matrix;
			USparseMatrix system_matrix;
			UVector pressure_w;
			UVector mass_vector;
			UVector f;

			UVector upwind_vector;

			Factorization<USparseMatrix, UVector> mass_matrix_inverse;
			double dt;
			double simulation_time;

			std::vector<double> box_min, box_max;
			std::vector<int> in_out_flow;
		};

		SteadyFractureFlowSimulation steady_flow_;
		USparseMatrix upwind_f_;

		Transport transport_m_;
		Transport transport_f_;
		
		bool preset_velocity_field_;

	};

}

#endif //UTOPIA_FRACTURE_FLOW_TRANSPORT_SIMULATION_HPP
