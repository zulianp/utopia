#ifndef UTOPIA_UI_HPP
#define UTOPIA_UI_HPP

#include <memory>
#include <string>

#include "utopia_InputStream.hpp"
#include "utopia_Path.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

	class Path;
	std::unique_ptr<InputStream> open_istream(const Path &path);
}

#endif //UTOPIA_UI_HPP
