#ifndef UTOPIA_OLD_CONFIGURATION_HPP
#define UTOPIA_OLD_CONFIGURATION_HPP


#include "utopia.hpp"
#include <array>

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"

#include "utopia_Input.hpp"
#include "utopia_AffineTransform.hpp"
#include "utopia_ProductFunctionSpace.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include "utopia_GaitCycle.hpp"

namespace utopia {
	class OldConfiguration final : public GaitCycle::Configuration {
	public:
		using Point2d = std::array<double, 2>;
		using Point3d = std::array<double, 3>;
		using Fun2d = std::function<std::array<double, 2>(const std::array<double, 2> &p)>;
		using Fun3d = std::function<std::array<double, 3>(const std::array<double, 3> &p)>;

		OldConfiguration();
		void read(Input &is) override;
		void update(const int time_step) override;
		void displacement_and_forces(
			ProductFunctionSpace<LibMeshFunctionSpace> &space,
			UVector &displacement,
			UVector &forces) const override;

		void init(ProductFunctionSpace<LibMeshFunctionSpace> &space) override;

		inline int n_steps() const override
		{
			return n_time_steps;
		}

		void describe(std::ostream &os) const override
		{
			os << "n_time_steps: " << n_time_steps << std::endl;
			std::cout << dt_ << std::endl;
			std::cout << t << std::endl;
			std::cout << t_end << std::endl;
		}

		double dt() const override;

		inline static std::string name()
		{
			return "old";
		}

		class Rotation {
		public:
			int block;
			int n_dims;
			char axis;
			double begin_angle_degree;
			double end_angle_degree;
			double d_angle;

			AffineTransform trafo;

			int from_step;
			int to_step;

			bool active;

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

			int from_step;
			int to_step;

			bool active;

			AffineTransform trafo;

			Translation();

			void init(const int n_dims, const int n_steps);
			void update(const int step, const double t);
		};


		int n_time_steps;
		double dt_;
		double t;
		double t_end;

		std::vector<Rotation> rotations;
		std::vector<Translation> translations;

	};
}

#endif //UTOPIA_OLD_CONFIGURATION_HPP
