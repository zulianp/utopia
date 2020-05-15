#ifndef UTOPIA_CL_ERROR_HPP
#define UTOPIA_CL_ERROR_HPP

#include <iostream>
#include "cl.hpp"

namespace utopia {
    bool check_cl_error(const cl_int code);
}

#endif  // UTOPIA_CL_ERROR_HPP
