#ifndef UTOPIA_ENUMS_HPP
#define UTOPIA_ENUMS_HPP

namespace utopia {

	enum WriteMode {
	    AUTO  = 0, //unsafe for the moment depends on backend implementation
	    LOCAL = 1,
	    GLOBAL_INSERT = 2,
	    GLOBAL_ADD    = 3
	};
	
}

#endif //UTOPIA_ENUMS_HPP
