# ParMetis_FOUND       - Do not attempt to use if "no" or undefined
# ParMetis_INCLUDES    - the ParMETIS includes
# ParMetis_LIBRARIES   - List of fully qualified libraries to link against

# if(NOT ParMetis_FOUND)
    find_path(ParMetis_INCLUDES parmetis.h
        PATHS
        ${ParMetis_INCLUDES}
        ${ParMetis_INCLUDE_PATH}
        ${ParMetis_ROOT}/include
        ${ParMetis_DIR}/include
        $ENV{PARMETIS_INCLUDE_PATH}
        $ENV{PARMETIS_ROOT}/include
        $ENV{PARMETIS_DIR}/include
    	/usr/local/include
    	/usr/include
    	/usr/include/metis
    )

    find_library(ParMetis_LIBRARY parmetis
        PATHS
        ${ParMetis_LIBRARY}
        ${ParMetis_LIBRARY_PATH}
        ${ParMetis_DIR}
        $ENV{PARMETIS_LIBRARY_PATH}
        $ENV{PARMETIS_ROOT}
        $ENV{PARMETIS_DIR}
    	/usr/local
    	/usr	
        PATH_SUFFIXES lib
    )

    message(STATUS "ParMetis_INCLUDES=${ParMetis_INCLUDES}\nParMetis_LIBRARY=${ParMetis_LIBRARY}")
    	
    if ( ParMetis_INCLUDES )	
    	if( ParMetis_LIBRARY )
    		set( ParMetis_LIBRARIES ${ParMetis_LIBRARY})
    		set( ParMetis_FOUND TRUE )
    	endif ()
    endif ()


    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(ParMetis DEFAULT_MSG
                                      ParMetis_INCLUDES ParMetis_LIBRARIES)

    mark_as_advanced ( ParMetis_FOUND ParMetis_INCLUDES ParMetis_LIBRARIES )
# endif()