# TODO see if we can use Trilinos' target to retrieve Intrepid2

find_package(Intrepid2 REQUIRED)

utopiafe_add_module(utopia_intrepid2 ${CMAKE_CURRENT_SOURCE_DIR} ".")

target_link_directories(utopia_intrepid2 PUBLIC ${Intrepid2_LIBRARY_DIRS})
target_link_libraries(utopia_intrepid2 PUBLIC utopia_fe_core Utopia::utopia)
target_link_libraries(utopia_intrepid2 PUBLIC ${Intrepid2_LIBRARIES})

target_include_directories(utopia_intrepid2
                           PUBLIC ${Intrepid2_TPL_INCLUDE_DIRS})

target_include_directories(utopia_intrepid2 PUBLIC ${Intrepid2_INCLUDE_DIRS})

set(UTOPIA_WITH_INTREPID2 TRUE)

set(UTOPIA_WITH_INTREPID2
    TRUE
    PARENT_SCOPE)
