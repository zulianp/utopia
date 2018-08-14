#ifndef UTOPIA_WEAR_HPP
#define UTOPIA_WEAR_HPP

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"

#include "utopia_FEForwardDeclarations.hpp"

#include "utopia.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_Mechanics.hpp"


#include <vector>
#include <cmath>

namespace utopia {

	void apply_displacement(
		const DVectord &displacement_increment,
		const libMesh::DofMap &dof_map,
		libMesh::MeshBase &mesh);

	class Wear {
	public:

		inline void update(
			const double dt,
			const DVectord &sliding_distance,
			const DVectord &normal_stress)
		{
			if(empty(wear)) {
				wear = local_zeros(local_size(normal_stress));
			}

			wear += (dt * wear_coefficient) * abs(e_mul(sliding_distance, normal_stress));
		}

		void compute_displacement(
			ProductFunctionSpace<LibMeshFunctionSpace> &V,
			const std::vector<int> &boundary_tags,
			DVectord &wear_induced_displacement
			);

		void modify_geometry(
			ProductFunctionSpace<LibMeshFunctionSpace> &V,
			const std::vector<int> &boundary_tags);

		void init_aux_system(
			libMesh::EquationSystems &es,
			libMesh::Order order = libMesh::FIRST);

		void update_aux_system(
			const int main_system_number,
			const MechanicsState &state,
			const Contact &contact,
			const double dt,
			libMesh::EquationSystems &es);

	private:
		DVectord wear;
		double wear_coefficient;
		double extrapolation_factor;

		std::vector<int> var_num_aux;
		unsigned int param_sys_number;

		std::vector<double> total_wear;

		//buffers
		DVectord wear_induced_displacement;
		DVectord is_normal_component;
		DVectord normals;
		DSMatrixd trafo;
	};
}

#endif //UTOPIA_WEAR_HPP
