# PARMETIS_FOUND       - Do not attempt to use if "no" or undefined
# PARMETIS_INCLUDES    - the ParMETIS includes
# PARMETIS_LIBRARIES   - List of fully qualified libraries to link against
	
if ( NOT DEFINED PARMETIS_INCLUDES )
    find_path(PARMETIS_INCLUDES parmetis.h
        PATHS
        ${PARMETIS_INCLUDE_PATH}
        ${PARMETIS_ROOT}/include
        $ENV{PARMETIS_INCLUDE_PATH}
        $ENV{PARMETIS_ROOT}/include
    	/usr/local/include
    	/usr/include
    	/usr/include/metis
    )
endif ()

find_library(PARMETIS_LIBRARY parmetis
    PATHS
    ${PARMETIS_LIBRARY_PATH}
    $ENV{PARMETIS_LIBRARY_PATH}
    $ENV{PARMETIS_ROOT}
	/usr/local
	/usr	
    PATH_SUFFIXES lib
)
	
if ( PARMETIS_INCLUDES )	
	if( PARMETIS_LIBRARY )
		set( PARMETIS_LIBRARIES ${PARMETIS_LIBRARY})
		set( PARMETIS_FOUND TRUE )
        
	endif ()
endif ()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARMETIS DEFAULT_MSG
                                  PARMETIS_INCLUDES PARMETIS_LIBRARIES)


mark_as_advanced ( PARMETIS_FOUND PARMETIS_INCLUDES PARMETIS_LIBRARIES )
