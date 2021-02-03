utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")

# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
if(UTOPIA_ENABLE_BOOST)
    find_package(Boost)
    find_package(BoostExtended)
    if(Boost_FOUND)
        target_link_libraries(utopia_fe PUBLIC Boost::boost)
        set(UTOPIA_WITH_BOOST
            TRUE
            PARENT_SCOPE)
    else()
        message(WARNING "Boost not found.")
        set(UTOPIA_WITH_BOOST
            FALSE
            PARENT_SCOPE)
    endif()
else()
    set(UTOPIA_WITH_BOOST
        FALSE
        PARENT_SCOPE)
endif()
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
