#include "utopia_ui.hpp"
#include "utopia_Path.hpp"
#include "utopia_XMLInput.hpp"
#include "utopia_JSONInput.hpp"
#include "utopia_Convertible.hpp"
#include "utopia_InputParameters.hpp"


namespace utopia {
	template class Convertible<double>;
	template class Convertible<float>;
	template class Convertible<long>;
	template class Convertible<int>;
	template class Convertible<bool>;
	template class Convertible<std::string>;


    std::unique_ptr<Input> open_istream(const Path &path)
	{
		if(path.extension() == "xml") {
			// auto ret = make_unique<XMLInput>();
			// ret->open(path);
			// return ret;
			
			//portable version when compuling with nvcc
			auto ret = new XMLInput();
			auto ret_ptr = std::unique_ptr<Input>(ret);
			ret->open(path);
			return ret_ptr;
			
		} else if(path.extension() == "json") {

			auto ret = new JSONInput();
			auto ret_ptr = std::unique_ptr<Input>(ret);
			ret->open(path);
			return ret_ptr;
		} else {
			std::cerr << "[Error] format not supported" << std::endl;
			return nullptr;
		}
	}
}
