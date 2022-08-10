find_package(ParMetis)

if(NOT PARMETIS_FOUND)
    message(FATAL_ERROR "Unable to find ParMetis")
    # include(FetchContent)
    # message(STATUS "Fetching parmetis, since it could not be found.")

    # if(UTOPIA_DEPENDENCIES_DIR)
    #     set(PARMETIS_INSTALL_DIR ${UTOPIA_DEPENDENCIES_DIR}/PARMETIS)
    # else()
    #     set(PARMETIS_INSTALL_DIR ${CMAKE_SOURCE_DIR}/../external/PARMETIS)
    # endif()

    # set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage")
    # set(PARMETIS_URL "https://github.com/KarypisLab/ParMETIS")

    # FetchContent_Declare(
    #     parmetis
    #     GIT_REPOSITORY ${PARMETIS_URL}
    #     # GIT_TAG origin/sampler
    #     GIT_TAG origin/development
    # )
    
    # FetchContent_MakeAvailable(parmetis)

    # add_library(ParMetis::parmetis ALIAS parmetis)
    # target_link_libraries(utopia PUBLIC parmetis)

else()
    target_include_directories(utopia PUBLIC ${PARMETIS_INCLUDES})
    target_link_libraries(utopia PUBLIC ${PARMETIS_LIBRARIES})

    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${PARMETIS_INCLUDES})
    list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${PARMETIS_LIBRARIES})

    set(UTOPIA_THIRDPARTY_INCLUDES
        ${UTOPIA_THIRDPARTY_INCLUDES}
        PARENT_SCOPE)

    set(UTOPIA_THIRDPARTY_LIBRARIES
        ${UTOPIA_THIRDPARTY_LIBRARIES}
        PARENT_SCOPE)
endif()


set(UTOPIA_WITH_PARMETIS TRUE)
set(UTOPIA_WITH_PARMETIS TRUE PARENT_SCOPE)

list(
    APPEND
    PARMETIS_MODULES
    .
    )

utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${PARMETIS_MODULES}")
