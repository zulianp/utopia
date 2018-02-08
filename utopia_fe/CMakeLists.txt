cmake_minimum_required(VERSION 2.8)
project(utopia_fe)

#To kill the policy warning  (maybe not a good idea yet)
set(CMAKE_MACOSX_RPATH 1)

option(ENABLE_SANITIZER "check for memory access problems" OFF)
option(ENABLE_BOOST "enable boost features" OFF)
option(ENABLE_MOONOLITH_PROFILER "enable moonolith profiling capabilities" OFF)

#################################################################################
#################################################################################
#################################################################################

set(UTOPIA_FE_ROOT_PATH ${CMAKE_SOURCE_DIR})
set(UTOPIA_FE_CMAKES_PATH ${UTOPIA_FE_ROOT_PATH}/libcmake)


if (ENABLE_SANITIZER AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")
		SET(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls ")
	endif()
endif()

################################## Include cmake functions and find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_FE_ROOT_PATH}/cmake")
file(GLOB_RECURSE FUNC_LIB "{UTOPIA_FE_ROOT_PATH}/cmake/utils" )
foreach(funcLib ${FUNC_LIB})
	include(${funcLib})
endforeach(funcLib)

find_package(Utopia)
if(UTOPIA_FOUND)
    list(APPEND CMAKE_CXX_FLAGS "-std=c++11")
    include_directories(${UTOPIA_INCLUDES})
    add_definitions(${UTOPIA_DEFS})
else()
    message(FATAL_ERROR "Utopia not found.")
endif()

include_directories(rapidxml)


########################################################################
function(find_project_files rootPath dirPaths headers sources)
	SET(verbose TRUE)

	SET(theaders ${${headers}})
	SET(tsources ${${sources}})

	SET(ROOT_PATH ${${rootPath}})
	SET(DIR_PATHS ${${dirPaths}})

	FOREACH(INCLUDE_PATH ${DIR_PATHS})
		INCLUDE_DIRECTORIES(${ROOT_PATH}/${INCLUDE_PATH})

		file(GLOB TEMP_HPPSRC "${ROOT_PATH}/${INCLUDE_PATH}/*.cpp" )   
		file(GLOB TEMP_SRC "${ROOT_PATH}/${INCLUDE_PATH}/*.c" )    
		file(GLOB TEMP_HPPDR "${ROOT_PATH}/${INCLUDE_PATH}/*.hpp" )
		file(GLOB TEMP_HDR "${ROOT_PATH}/${INCLUDE_PATH}/*.h" )
		
		SOURCE_GROUP(${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR}; ${TEMP_HPPSRC}; ${TEMP_SRC}; ${TEMP_UI})
		
		SET(tsources ${tsources}; ${TEMP_SRC}; ${TEMP_HPPSRC})
		SET(theaders ${theaders}; ${TEMP_HDR}; ${TEMP_HPPDR}) 
	ENDFOREACH(INCLUDE_PATH)  

	SET(${headers} ${theaders} PARENT_SCOPE)
	SET(${sources} ${tsources} PARENT_SCOPE)
endfunction()
########################################################################


list(APPEND UTOPIA_FE_MODULES
	tests
	transfer
	visualization
	backend/libmesh
	backend/homemade
	intersection
	materials
	core
	rapidxml
 )

FOREACH(MODULE ${UTOPIA_FE_MODULES})
		add_subdirectory("${MODULE}")
ENDFOREACH(MODULE)  


include_directories(${CMAKE_BINARY_DIR})
configure_file(${UTOPIA_FE_ROOT_PATH}/utopia_fe_config.hpp.in ${CMAKE_BINARY_DIR}/utopia_fe_config.hpp)
install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.hpp DESTINATION include)



#create main library
list(APPEND TOP_HEADERS utopia_fe.hpp)
list(APPEND TOP_SOURCES utopia_fe.cpp)


add_library(utopia_fe STATIC ${TOP_SOURCES} ${TOP_HEADERS})
target_link_libraries(utopia_fe ${UTOPIA_LIBRARIES} ${UTOPIA_FE_LIBRARIES} ${TRANSFER_LIB})

add_executable(utopia_fe_exec  main.cpp ${TOP_HEADERS})
# target_link_libraries(utopia_fe_exec utopia_fe ${UTOPIA_FE_LIBRARIES})
# these next three lines fix some issues on linux
target_link_libraries(utopia_fe_exec utopia_fe)
target_link_libraries(utopia_fe_exec ${UTOPIA_FE_LIBRARIES})
target_link_libraries(utopia_fe_exec ${UTOPIA_LIBRARIES})


install(TARGETS utopia_fe
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(TARGETS utopia_fe_exec
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${TOP_HEADERS} DESTINATION include)




#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

function(append_list_to_string_for_makefile OUTPUT_STRING PRE_LIST POST_LIST ITEM_PREFIX LIST_OF_STRINGS)
	#http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro

	LIST(REMOVE_DUPLICATES LIST_OF_STRINGS)

	SET(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

	FOREACH(LIST_EL ${LIST_OF_STRINGS}) 
		
		string(REGEX MATCH "(.*)framework" IS_FRAME_WORK "${LIST_EL}")
		if(IS_FRAME_WORK)
			get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
			string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH "${LIST_EL}")
			string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME "${FRAMEWORK_NAME}")

			SET(TEMPSTR "${TEMPSTR}-F${FRAMEWORK_PATH} -framework ${FRAMEWORK_NAME} ")
		else()
	    	SET(TEMPSTR "${TEMPSTR}${ITEM_PREFIX}${LIST_EL} ")
	    endif()
	ENDFOREACH(LIST_EL) 

	SET(TEMPSTR "${TEMPSTR} ${POST_LIST}")
	SET(${OUTPUT_STRING} ${TEMPSTR} PARENT_SCOPE)
endfunction()

function(append_list_to_string_for_makefile_with_postfix OUTPUT_STRING PRE_LIST POST_LIST ITEM_PREFIX ITEM_POSTFIX LIST_OF_STRINGS)
	#http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro

	LIST(REMOVE_DUPLICATES LIST_OF_STRINGS)

	SET(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

	FOREACH(LIST_EL ${LIST_OF_STRINGS}) 
		
		string(REGEX MATCH "(.*)framework" IS_FRAME_WORK "${LIST_EL}")
		if(IS_FRAME_WORK)
			get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
			string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH "${LIST_EL}")
			string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME "${FRAMEWORK_NAME}")

			SET(TEMPSTR "${TEMPSTR}-F${FRAMEWORK_PATH} -framework ${FRAMEWORK_NAME} ")
		else()
	    	SET(TEMPSTR "${TEMPSTR}${ITEM_PREFIX}${LIST_EL}${ITEM_POSTFIX} ")
	    endif()
	ENDFOREACH(LIST_EL) 

	SET(TEMPSTR "${TEMPSTR} ${POST_LIST}")
	SET(${OUTPUT_STRING} ${TEMPSTR} PARENT_SCOPE)
endfunction()

#################################################################################

function(append_list_to_string_for_cmake OUTPUT_STRING PRE_LIST POST_LIST LIST_OF_STRINGS)
#http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro
	LIST(REMOVE_DUPLICATES LIST_OF_STRINGS)

	SET(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

	FOREACH(LIST_EL ${LIST_OF_STRINGS}) 
	    SET(TEMPSTR "${TEMPSTR}\t\"${LIST_EL}\"\n")
	ENDFOREACH(LIST_EL) 

	SET(TEMPSTR "${TEMPSTR} ${POST_LIST}")
	SET(${OUTPUT_STRING} ${TEMPSTR} PARENT_SCOPE)
endfunction()



# Generate files useful for dependent build systems
set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n#Includes the content of utopia_config.makefile generated by utopia\n\n")

list(APPEND UTOPIA_FE_INCLUDES ${UTOPIA_INCLUDES})
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_FE_INCLUDES=-I${CMAKE_INSTALL_PREFIX}/include -I${LIBMESH_INCLUDE_DIR} " " " "-I" "${UTOPIA_FE_INCLUDES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" "\n" "-I" "${MOONOLITH_INCLUDES}")

append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_FE_LIBRARIES=-L${CMAKE_INSTALL_PREFIX}/lib -lutopia_fe ${LIBMESH_LIBRARIES} " " " "-l" "${UTOPIA_FE_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " "" "${UTOPIA_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " "" "${MOONOLITH_LIBRARIES}")

set(CONFIG_CONTENT "${CONFIG_CONTENT}\n")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_FE_DEFS=" "\n" "" "${UTOPIA_DEFS}")


set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_FE_CXX_FLAGS=${UTOPIA_CXX_FLAGS}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_FE_CXX_COMPILER=${UTOPIA_CXX_COMPILER}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_FE_C_COMPILER=${UTOPIA_C_COMPILER}\n")

file(WRITE ${CMAKE_BINARY_DIR}/utopia_fe_config.makefile ${CONFIG_CONTENT})

set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_INCLUDES\n\t${CMAKE_INSTALL_PREFIX}/include\n" 		"" "${UTOPIA_FE_INCLUDES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "" "" "${LIBMESH_INCLUDE_DIR}")
append_list_to_string_for_cmake(CONFIG_CONTENT "" ")\n\n"  ${MOONOLITH_INCLUDES})
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_LIBRARIES\n\t-L${CMAKE_INSTALL_PREFIX}/lib" "" "\n\t\"-l" "\"" "${UTOPIA_FE_LIBRARIES}")
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "\n\t-L${CMAKE_INSTALL_PREFIX}/lib" "" "\n\t\"" "\"" "${MOONOLITH_LIBRARIES}")
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "" "" "\n\t\"" "\"" "${LIBMESH_LIBRARIES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "\n" "\n)\n\n" "${UTOPIA_LIBRARIES}")

append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_DEFS\n" 			"" "${UTOPIA_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "" 			"\n)\n\n" "${UTOPIA_FE_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_CXX_FLAGS\n"  	"\n)\n\n" "${CMAKE_CXX_FLAGS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_C_COMPILER\n"  	"\n)\n\n" "${CMAKE_C_COMPILER}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_CXX_COMPILER\n" 	"\n)\n\n" "${CMAKE_CXX_COMPILER}")


file(WRITE ${CMAKE_BINARY_DIR}/utopia_fe_config.cmake ${CONFIG_CONTENT})

install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.makefile DESTINATION config)
install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.cmake 	 DESTINATION config)

