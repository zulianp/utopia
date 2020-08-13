find_package(Intrepid2)
if(INTREPID2_FOUND)
    utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")
    target_link_libraries(utopia_fe PUBLIC intrepid2)
    set(WITH_INTREPID2 TRUE PARENT_SCOPE)
else()
    message(WARNING "could not find Intrepid2")
endif()

