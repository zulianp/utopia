#include "utopia_Input.hpp"
#include "utopia_ui.hpp"
#include <exception>

namespace utopia {
	bool Configurable::import(const Path &path)
	{
		try {
			auto istr = open_istream(path.to_string());
			if(!istr) {
				return false;
			}

			read(*istr);
			return true;
		} catch(const std::exception &ex) {
			std::cerr << "[Error] " << ex.what() << std::endl;
			return false;
		}
		
	}

	bool Configurable::import(
		const std::string &key,
		const Path &path)
	{
		try {
			auto istr = open_istream(path.to_string());
			if(!istr) {
				return false;
			}

			istr->get(key, *this);
			return true;
		} catch(const std::exception &ex) {
			std::cerr << "[Error] " << ex.what() << std::endl;
			return false;
		}
	}

	void Configurable::print_usage(std::ostream &) const {}
}