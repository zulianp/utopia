# interoperability

list(APPEND INTEROP_MODULES .)

set(INTEROP_DEP "")
set(INTEROP_EXISTS FALSE)

# message(STATUS "${UTOPIA_ENABLE_MOONOLITH} == ${UTOPIA_ENABLE_STK}")

if(UTOPIA_ENABLE_MOONOLITH)
    if(UTOPIA_ENABLE_STK)
        list(APPEND INTEROP_DEP utopia_moonolith utopia_stk)
        list(APPEND INTEROP_MODULES moonolith_stk)
        set(INTEROP_EXISTS TRUE)
    endif()

    if(UTOPIA_ENABLE_LIBMESH)
        list(APPEND INTEROP_DEP utopia_moonolith utopia_fe)
        list(APPEND INTEROP_MODULES moonolith_libmesh)
        set(INTEROP_EXISTS TRUE)
    endif()

endif()

if(UTOPIA_ENABLE_STK)
    if(UTOPIA_ENABLE_INTREPID2)
        list(APPEND INTEROP_DEP utopia_intrepid2 utopia_stk)
        list(APPEND INTEROP_MODULES stk_intrepid2)
        set(INTEROP_EXISTS TRUE)
    endif()
endif()

if(INTEROP_EXISTS)
    utopiafe_add_module(utopia_interop ${CMAKE_CURRENT_SOURCE_DIR}
                        "${INTEROP_MODULES}")

    target_link_libraries(utopia_interop PUBLIC Utopia::utopia)
    target_link_libraries(utopia_interop PUBLIC ${INTEROP_DEP})

    target_include_directories(utopia_interop BEFORE
                               PUBLIC ${CMAKE_SOURCE_DIR}/base)

    target_include_directories(utopia_interop BEFORE
                               PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    target_include_directories(utopia_interop BEFORE PUBLIC ${CMAKE_BINARY_DIR})

    set(INTEROP_LIBRARIES
        ${INTEROP_LIBRARIES}
        PARENT_SCOPE)
    set(INTEROP_INCLUDES
        ${INTEROP_INCLUDES}
        PARENT_SCOPE)

endif()
