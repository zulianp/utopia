#include "utopia_CLError.hpp"
#include <iostream>

#define UTOPIA_CL_CHECK_CASE(ERROR_CODE)    \
    case ERROR_CODE: {                            \
        message = #ERROR_CODE;                    \
        std::cerr << "[Error] " << message << std::endl;        \
        return false;                            \
    }

namespace utopia {
    bool check_cl_error(const cl_int code) {
        std::string message;

        switch (code) {
            case CL_SUCCESS: {
                message = "";
                return true;
            }

            UTOPIA_CL_CHECK_CASE(CL_INVALID_PLATFORM)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_DEVICE_TYPE)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_VALUE)
            UTOPIA_CL_CHECK_CASE(CL_DEVICE_NOT_FOUND)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_COMMAND_QUEUE)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_OPERATION)
            UTOPIA_CL_CHECK_CASE(CL_COMPILER_NOT_AVAILABLE)
            UTOPIA_CL_CHECK_CASE(CL_OUT_OF_HOST_MEMORY)
            UTOPIA_CL_CHECK_CASE(CL_BUILD_PROGRAM_FAILURE)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_DEVICE)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_HOST_PTR);

                //clSetKernelArg
            UTOPIA_CL_CHECK_CASE(CL_INVALID_KERNEL)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_ARG_INDEX)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_ARG_VALUE)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_MEM_OBJECT)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_ARG_SIZE)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_SAMPLER)


            UTOPIA_CL_CHECK_CASE(CL_MAP_FAILURE);
            UTOPIA_CL_CHECK_CASE(CL_INVALID_GLOBAL_WORK_SIZE);

            UTOPIA_CL_CHECK_CASE(CL_INVALID_GL_OBJECT)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_CONTEXT)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_QUEUE_PROPERTIES)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_KERNEL_NAME)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_KERNEL_DEFINITION)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_EVENT_WAIT_LIST)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_WORK_DIMENSION)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_WORK_GROUP_SIZE)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_WORK_ITEM_SIZE)
            UTOPIA_CL_CHECK_CASE(CL_MEM_OBJECT_ALLOCATION_FAILURE)
            UTOPIA_CL_CHECK_CASE(CL_INVALID_KERNEL_ARGS)

            default: {
                std::cerr << "unknown return value for " << code << std::endl;
                return false;
            }
        }

        return false;
    }
}

#undef UTOPIA_CL_CHECK_CASE
