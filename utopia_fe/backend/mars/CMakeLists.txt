if(UTOPIA_ENABLE_MARS)
    find_package(MARS REQUIRED)

    message(STATUS "MARS_LIBRARIES=${MARS_LIBRARIES}")

    utopiafe_add_module(utopia_mars ${CMAKE_CURRENT_SOURCE_DIR} ".;assemblers;io;fe")
    target_include_directories(utopia_mars PUBLIC ${MARS_INCLUDES})
    target_link_libraries(utopia_mars INTERFACE ${MARS_LIBRARIES})
    target_link_libraries(utopia_mars PUBLIC Utopia::utopia utopia_fe_core utopia_fe_kokkos)

    # target_compile_features(utopia_mars PUBLIC cxx_std_17)

    set(UTOPIA_WITH_MARS
        TRUE
        PARENT_SCOPE)

    if(UTOPIA_ENABLE_MARS_VTK)
        find_package(
            VTK
            COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersGeneral vtkIOXML vtkIOParallel vtkIOParallelXML
            REQUIRED)

        if(VTK_FOUND)
            target_link_libraries(utopia_mars PUBLIC ${VTK_LIBRARIES})

            set(UTOPIA_WITH_VTK
                TRUE
                PARENT_SCOPE)

        else()
            message(FATAL_ERROR "VTK not FOUND")
        endif()
    endif()

endif()
