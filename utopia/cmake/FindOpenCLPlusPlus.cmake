FIND_PACKAGE(OpenCL)
IF(OPENCL_FOUND)
    FIND_PATH(OPEN_CL_PLUS_PLUS_SRC_PATH opencl/cl.hpp
    HINT ../external
    DOC "The OpenCL c++ wrapper"
    )

    # IF(OPEN_CL_PLUS_PLUS_SRC_PATH)
    # 	SET(OPENCL_INCLUDE_DIRS; ${OPENCL_INCLUDE_DIRS} ${OPEN_CL_PLUS_PLUS_SRC_PATH}/opencl)
    # 	ELSE()
    # 		MESSAGE(WARNING "OpenCL++ wrapper not available")
    # ENDIF()
ENDIF()