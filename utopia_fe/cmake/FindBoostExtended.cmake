FIND_PACKAGE(Boost COMPONENTS thread system)

IF(NOT Boost_FOUND)
	IF(WIN32)

		FIND_PATH(BOOST_ROOT_PATH  "boost/thread.hpp"
				  HINTS $ENV{BOOST_ROOT}
				  		"C:/boost"
				  		"C:/local/boost_1_54_0"
		)
		
		

		IF(BOOST_ROOT_PATH)
			SET(Boost_INCLUDE_DIRS ${BOOST_ROOT_PATH})


			FIND_PATH(Boost_LINK_DIRECTORY "libboost_thread-$ENV{VISUAL_STUDIO_COMPILER}-mt-$ENV{BOOST_VERSION_NUMBER}.lib"
					  HINTS "$ENV{BOOST_LIB_DIR}"
					  		"${BOOST_ROOT_PATH}/lib32-msvc-$ENV{VISUAL_STUDIO_COMPILER}"
					  		"${BOOST_ROOT_PATH}/stage/lib"
					  		"${BOOST_ROOT_PATH}/lib32-msvc-11.0"
					  		"${BOOST_ROOT_PATH}/lib32-msvc-12.0"
			)
			
			IF(Boost_LINK_DIRECTORY)
				SET(Boost_LIBRARY_DIRS ${Boost_LINK_DIRECTORY})

				IF(DEBUG)
					SET(Boost_LIBRARIES ${Boost_LIBRARIES}; libboost_thread-$ENV{VISUAL_STUDIO_COMPILER}-mt-gd-$ENV{BOOST_VERSION_NUMBER}; libboost_system-$ENV{VISUAL_STUDIO_COMPILER}-mt-gd-$ENV{BOOST_VERSION_NUMBER}; boost_date_time-$ENV{VISUAL_STUDIO_COMPILER}-mt-gd-$ENV{BOOST_VERSION_NUMBER})
					# SET(Boost_LIBRARIES ${Boost_LIBRARIES}; libboost_thread-vc110-mt-gd-1_54; libboost_system-vc110-mt-gd-1_54; boost_date_time-vc110-mt-gd-1_54)
				ELSE()
					SET(Boost_LIBRARIES ${Boost_LIBRARIES}; libboost_thread-$ENV{VISUAL_STUDIO_COMPILER}-mt-$ENV{BOOST_VERSION_NUMBER}; libboost_system-$ENV{VISUAL_STUDIO_COMPILER}-mt-$ENV{BOOST_VERSION_NUMBER}; boost_date_time-$ENV{VISUAL_STUDIO_COMPILER}-mt-$ENV{BOOST_VERSION_NUMBER})
					# SET(Boost_LIBRARIES ${Boost_LIBRARIES}; libboost_thread-vc110-mt-1_54; libboost_system-vc110-mt-1_54; boost_date_time-vc110-mt-1_54)
				ENDIF()

				MESSAGE(WARNING "[Warning] Using some default hardcoded libs")
				SET(Boost_FOUND TRUE)
			ENDIF()
		ENDIF()
	ENDIF()
ENDIF()