# find dependencies
find_package(Trilinos)

if(TRILINOS_FOUND)
    ##### These lines are needed so the utopia-config.makefile will be correctly built. The UtopiaTargets.cmake file is populated automatically by cmake #########
   foreach(LIB ${Trilinos_LIBRARY_DIRS})
       list(APPEND UTOPIA_THIRDPARTY_LIBRARIES -L${LIB})
   endforeach(LIB)

   list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${Trilinos_INCLUDE_DIRS})
   # list(APPEND UTOPIA_LIBRARIES ${Trilinos_LIBRARY_DIRS})
   list(APPEND Trilinos_all_libs ${Trilinos_LIBRARIES})
   list(APPEND Trilinos_all_libs ${Trilinos_TPL_LIBRARIES})
   list(REVERSE Trilinos_all_libs)
   list(REMOVE_DUPLICATES Trilinos_all_libs)
   list(REVERSE Trilinos_all_libs)
   foreach(LIB ${Trilinos_all_libs})
       if(EXISTS ${LIB})
           list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${LIB})
       else()
           list(APPEND UTOPIA_THIRDPARTY_LIBRARIES -l${LIB})
       endif()
   endforeach()
   list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${Trilinos_TPL_INCLUDE_DIRS})
   ############### End of population of variables for utopia-config.makefile #########################3

    set(UTOPIA_TRILINOS_DEPS Kokkos::kokkos)
    # Unfortunately trilinos libraries are not namespaced (yet), we use the bare target names
    list(APPEND Utopia_Trilinos_possible_packages Amesos2 Belos Ifpack2 MueLu)
    foreach(package ${Utopia_Trilinos_possible_packages})
        # e.g Amesos2: WITH_TRILINOS_AMESOS2=TRUE if package is found in Trilinos_PACKAGE_LIST
        string(TOUPPER ${package} packageUpper)
        string(TOLOWER ${package} packageLower)
        list(FIND Trilinos_PACKAGE_LIST ${package} PACKAGE_FOUND)
        if (NOT PACKAGE_FOUND EQUAL -1)
            set(WITH_TRILINOS_${packageUpper} TRUE PARENT_SCOPE)
            if(${package} STREQUAL "MueLu")
                # ugly hack, but we need to link with muelu-adapters also
                # I cannot wait until Trilinos finally supports cmake targets correctly
                list(APPEND UTOPIA_TRILINOS_DEPS "muelu-adapters")
            endif()
            list(APPEND UTOPIA_TRILINOS_DEPS ${packageLower})
        endif()
    endforeach()

    set(WITH_TRILINOS TRUE)
    set(UTOPIA_TRILINOS TRUE PARENT_SCOPE)

    find_package(TpetraExt)
    if(TRILINOS_TPETRAEXT_FOUND)
        set(WITH_TRILINOS_TPETRAEXT TRUE PARENT_SCOPE)
    endif()
else()
    message(WARNING "[Warning] Trilinos not found")
    set(WITH_TRILINOS FALSE)
endif()

# set-up module
if(TRILINOS_FOUND)
    list(APPEND TRILINOS_MODULES
        .
        eval
        solvers
        solvers/linear
        solvers/nonlinear
        kokkos
        )
    utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${TRILINOS_MODULES}")

    # include manually thirtparty library include directories (I think this should be added with target_link_libraries(utopia_trilinos TrilinosPaackage), but it is not)
    target_include_directories(utopia PUBLIC ${Trilinos_TPL_INCLUDE_DIRS})
    target_link_libraries(utopia PUBLIC ${UTOPIA_TRILINOS_DEPS})

    # target_compile_definitions(utopia PUBLIC ${TRILINOS_DEFINITIONS})
    # target_compile_definitions(utopia PUBLIC ${Kokkos_CXX_FLAGS})

    set(UTOPIA_THIRDPARTY_INCLUDES ${UTOPIA_THIRDPARTY_INCLUDES} PARENT_SCOPE)
    set(UTOPIA_THIRDPARTY_LIBRARIES ${UTOPIA_THIRDPARTY_LIBRARIES} PARENT_SCOPE)
endif()
