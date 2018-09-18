#include "utopia_GaitCycle.hpp"
#include "utopia_AffineTransform.hpp"

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include "utopia_OldConfiguration.hpp"
#include "utopia_ISO14243_3.hpp"

namespace utopia {
	void GaitCycle::read(InputStream &is)
	{
		// std::string type = ISO14243_3::name();
		std::string type = OldConfiguration::name();
		is.read("type", type);

		auto it = conf_types_.find(type);

		if(it == conf_types_.end()) {
			std::cerr << "[Error] " << type << " not found. Using " << OldConfiguration::name() << std::endl;
			// conf_ = std::make_shared<ISO14243_3>();
			conf_ = std::make_shared<OldConfiguration>();
		} else {

			std::cout << "[Status] using type " << type << std::endl;
			conf_ = it->second;
		}

		conf_->read(is);
	}

	GaitCycle::GaitCycle()
	{
		conf_types_[OldConfiguration::name()] = std::make_shared<OldConfiguration>();
		conf_types_[ISO14243_3::name()]       = std::make_shared<ISO14243_3>();
	}

	void GaitCycle::update(const std::size_t time_step)
	{
		assert(conf_);
		conf_->update(time_step);
	}

	void GaitCycle::displacement_and_forces(
		ProductFunctionSpace<LibMeshFunctionSpace> &space,
		UVector &displacement,
		UVector &forces) const
	{
		conf_->displacement_and_forces(space, displacement, forces);
	}
}
