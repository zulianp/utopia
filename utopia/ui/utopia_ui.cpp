#include "utopia_ui.hpp"
#include "utopia_Path.hpp"
#include "utopia_XMLStream.hpp"


namespace utopia {

	std::unique_ptr<InputStream> open_istream(const Path &path)
	{
		if(path.extension() == "xml") {
			// auto ret = make_unique<XMLInputStream>();
			// ret->open(path);
			// return ret;
			
			//portable version when compuling with nvcc
			auto ret = new XMLInputStream();
			auto ret_ptr = std::unique_ptr<InputStream>(ret);
			ret->open(path);
			return ret_ptr;
		} else {
			std::cerr << "[Error] format not supported" << std::endl;
			return nullptr;
		}
	}
}
