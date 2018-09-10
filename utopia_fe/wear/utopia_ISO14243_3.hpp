#ifndef UTOPIA_ISO14243_3_HPP
#define UTOPIA_ISO14243_3_HPP

#include "utopia.hpp"
#include <array>

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"

#include "utopia_InputStream.hpp"
#include "utopia_AffineTransform.hpp"
#include "utopia_ProductFunctionSpace.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include "utopia_GaitCycle.hpp"

namespace utopia {

	class ISO14243_3 final : public GaitCycle::Configuration {
	public:
		using Point2d = std::array<double, 2>;
		using Point3d = std::array<double, 3>;
		using Fun2d = std::function<std::array<double, 2>(const std::array<double, 2> &p)>;
		using Fun3d = std::function<std::array<double, 3>(const std::array<double, 3> &p)>;

		ISO14243_3() {}

		void read(InputStream &is) override;
		void update(const int time_step) override;
		void displacement_and_forces(
			ProductFunctionSpace<LibMeshFunctionSpace> &space,
			DVectord &displacement,
			DVectord &forces) const override;

		inline int n_steps() const override
		{
			return file_.n_rows();
		}

		inline static std::string name()
		{
			return "ISO14243_3";
		}

		inline double dt() const override {
			return dt_;
		}

	private:
		char flexion_extension_angle_axis_;			
		char ap_motion_axis_;
		char tibial_rotation_axis_;
		int axial_force_axis_;
		
		//blocks and side-sets
		int femural_block_;
		int tibial_block_;
		int axial_force_side_;


		//////////////////// Valus from ISO ////////////////
		//%
		int percentage_of_time_cycle_;
		//degrees
		double flexion_extension_angle_;

		//Newton
		double axial_force_;


		//mm
		double ap_motion_;

		//degrees
		double tibial_int_ext_rotation_;


		double dt_;
		
		CSV file_;


		inline bool read(const Path &path)
		{
			return file_.read(path);
		}
	};
}

#endif //UTOPIA_ISO14243_3_HPP
