#ifndef UTOPIA_FE_BASE_HPP
#define UTOPIA_FE_BASE_HPP

#include "utopia_fe_config.hpp"
#include "utopia.hpp"

namespace utopia {
#ifdef USE_UTOPIA_TRILINOS
	using USMatrix = TSMatrixd;
	using UVector  = TVectord;
#else
	using USMatrix = DSMatrixd;
	using UVector  = DVectord;
#endif //USE_UTOPIA_TRILINOS
}

#endif //UTOPIA_FE_BASE_HPP
