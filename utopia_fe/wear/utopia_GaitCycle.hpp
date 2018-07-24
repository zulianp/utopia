#ifndef UTOPIA_GAIT_CYCLE_HPP
#define UTOPIA_GAIT_CYCLE_HPP

#include "utopia.hpp"
#include <array>

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"


namespace utopia {

	class InputStream;

	class GaitCycle {
	public:
		using Point2d = std::array<double, 2>;
		using Point3d = std::array<double, 3>;
		using Fun2d = std::function<std::array<double, 2>(const std::array<double, 2> &p)>;
		using Fun3d = std::function<std::array<double, 3>(const std::array<double, 3> &p)>;

		void init(
			const int dim,
			InputStream &is);

		GaitCycle();

		void set_time_step(const std::size_t time_step);

		void toggle_dir();

		void init();

		void override_displacement(
			const libMesh::MeshBase &mesh,
			const libMesh::DofMap &dof_map,
			const int block_id_rot,
			const int block_id_trasl,
			DVectord &displacement) const;
		
		std::size_t n_time_steps;
		double dt;
		double t;
		double t_end;
		double start_angle_degree;
		double start_angle_radian;
		double angle_degree;
		double angle_radian;
		double d_angle;
		Fun2d rotate2;
		Fun3d rotate3;
		Fun2d translate2_y;
		Fun3d translate3_z;
		Fun2d zero2;
		Fun3d zero3;
		Fun3d bc34;
		bool negative_dir;
	};
}

#endif //UTOPIA_GAIT_CYCLE_HPP
