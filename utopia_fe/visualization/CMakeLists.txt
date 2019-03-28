list(APPEND LOCAL_HEADERS
    utopia_Socket.hpp
    )

list(APPEND LOCAL_SOURCES
    utopia_Socket.cpp
    )
add_utopiafe_module(utopia_fe_visualization utopia_fe_libmesh_backend transfer)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
if(ENABLE_BOOST)
    find_package(Boost)
    find_package(BoostExtended)
    if(Boost_FOUND)
            target_include_directories(utopia_fe_visualization PRIVATE ${Boost_INCLUDE_DIRS})
            link_directories(${Boost_LIBRARY_DIRS})
            target_link_libraries(utopia_fe_visualization INTERFACE ${Boost_LIBRARIES})
            # add_definitions("-DWITH_BOOST")
            set(WITH_BOOST TRUE PARENT_SCOPE)
    else()
        MESSAGE(WARNING "Boost not found.")
        set(WITH_BOOST FALSE PARENT_SCOPE)
    endif()
else()
    set(WITH_BOOST FALSE PARENT_SCOPE)
endif()
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
