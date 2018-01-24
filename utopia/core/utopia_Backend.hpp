//
// Created by Patrick Zulian on 18/05/15.
//

#ifndef utopia_utopia_BACKEND_HPP
#define utopia_utopia_BACKEND_HPP

#include <vector>
#include <memory>
#include <iostream>

#include "utopia_Base.hpp"

namespace utopia {

    template<typename Scalar, int BackendType>
    class Backend {
    public:
    	Backend()
    	{
    		static_assert(BackendType < HOMEMADE, "No Backend implemented");
    	}
    };
}
#endif //utopia_utopia_BACKEND_HPP
