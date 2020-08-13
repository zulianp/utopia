utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
if(ENABLE_BOOST)
    find_package(Boost)
    find_package(BoostExtended)
    if(Boost_FOUND)
        target_link_libraries(utopia_fe PUBLIC Boost::boost)
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
