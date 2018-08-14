#ifndef UTOPIA_GAIT_CYCLE_HPP
#define UTOPIA_GAIT_CYCLE_HPP

#include "utopia.hpp"
#include <array>

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"

#include "utopia_InputStream.hpp"
#include "utopia_AffineTransform.hpp"

namespace utopia {

	class InputStream;

	class GaitCycle : public Serializable {
	public:
		using Point2d = std::array<double, 2>;
		using Point3d = std::array<double, 3>;
		using Fun2d = std::function<std::array<double, 2>(const std::array<double, 2> &p)>;
		using Fun3d = std::function<std::array<double, 3>(const std::array<double, 3> &p)>;

		class Rotation {
		public:
			int block;
			int n_dims;
			char axis;
			double begin_angle_degree;
			double end_angle_degree;
			double d_angle;

			AffineTransform trafo;

			Rotation();

			void init(const int n_dims, const int n_steps);
			void update(const int step, const double t);
		};

		class Translation {
		public:
			int block;
			int n_dims;
			char axis;
			double begin_offset;
			double end_offset;
			double d_offset;

			AffineTransform trafo;

			Translation();

			void init(const int n_dims, const int n_steps);
			void update(const int step, const double t);
		};


		void read(InputStream &is) override;

		GaitCycle();

		void set_time_step(const std::size_t time_step);

		void toggle_dir();

		void init(int n_dims);
		void update_trafos(const int step, const double t);
		void init_trafos(const int n_dims, const int n_steps);

		void override_displacement(
			const libMesh::MeshBase &mesh,
			const libMesh::DofMap &dof_map,
			DVectord &displacement) const;


		void describe(std::ostream &os) const
		{
			os << "n_time_steps: " << n_time_steps << std::endl;
			std::cout << dt << std::endl;
			std::cout << t << std::endl;
			std::cout << t_end << std::endl;
		}

		int n_time_steps;
		double dt;
		double t;
		double t_end;
		// double start_angle_degree;
		// double start_angle_radian;
		// double angle_degree;
		// double angle_radian;
		// double d_angle;
		// Fun2d rotate2;
		// Fun3d rotate3;
		// Fun2d translate2_y;
		// Fun3d translate3_z;
		// Fun2d zero2;
		// Fun3d zero3;
		// Fun3d bc34;
		// bool negative_dir;

		std::vector<Rotation> rotations;
		std::vector<Translation> translations;
	};
}

#endif //UTOPIA_GAIT_CYCLE_HPP
