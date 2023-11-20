# ##############################################################################
function(find_project_files rootPath dirPaths headers sources)
    set(verbose TRUE)

    set(theaders ${${headers}})
    set(tsources ${${sources}})

    foreach(INCLUDE_PATH ${dirPaths})
        file(GLOB TEMP_HPPSRC "${rootPath}/${INCLUDE_PATH}/*.cpp")
        file(GLOB TEMP_SRC "${rootPath}/${INCLUDE_PATH}/*.c")
        file(GLOB TEMP_HPPDR "${rootPath}/${INCLUDE_PATH}/*.hpp")
        file(GLOB TEMP_HDR "${rootPath}/${INCLUDE_PATH}/*.h")

        source_group(
            ${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR}; ${TEMP_HPPSRC};
                                  ${TEMP_SRC}; ${TEMP_UI})

        set(tsources ${tsources}; ${TEMP_SRC}; ${TEMP_HPPSRC})
        set(theaders ${theaders}; ${TEMP_HDR}; ${TEMP_HPPDR})
    endforeach(INCLUDE_PATH)

    set(${headers}
        ${theaders}
        PARENT_SCOPE)
    set(${sources}
        ${tsources}
        PARENT_SCOPE)
endfunction()
# ##############################################################################

function(utopiafe_add_library libraryRootDir subDirs)
    set(LOCAL_HEADERS "")
    set(LOCAL_SOURCES "")
    find_project_files(${libraryRootDir} "${subDirs}" LOCAL_HEADERS
                       LOCAL_SOURCES)
    target_sources(utopia_fe PRIVATE ${LOCAL_SOURCES})
    install(FILES ${LOCAL_HEADERS} DESTINATION include)
    foreach(MODULE ${subDirs})
        target_include_directories(
            utopia_fe BEFORE
            PUBLIC $<BUILD_INTERFACE:${libraryRootDir}/${MODULE}>)
    endforeach(MODULE)
endfunction()

function(utopiafe_add_module module_name libraryRootDir subDirs)
    set(LOCAL_HEADERS "")
    set(LOCAL_SOURCES "")

    add_library(${module_name})
    find_project_files(${libraryRootDir} "${subDirs}" LOCAL_HEADERS
                       LOCAL_SOURCES)
    target_sources(${module_name} PRIVATE ${LOCAL_SOURCES})
    install(FILES ${LOCAL_HEADERS} DESTINATION include)
    foreach(MODULE ${subDirs})
        target_include_directories(
            ${module_name} BEFORE
            PUBLIC $<BUILD_INTERFACE:${libraryRootDir}/${MODULE}>)
    endforeach(MODULE)

    target_include_directories(
        ${module_name} BEFORE
        PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/base>)
    target_include_directories(
        ${module_name} BEFORE
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
    target_include_directories(${module_name} BEFORE
                               PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>)

    install(
        TARGETS ${module_name}
        EXPORT UtopiaFETargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

    target_link_libraries(utopia_fe PUBLIC ${module_name})
    list(APPEND UTOPIA_FE_LIBRARIES "${module_name}")
    set(UTOPIA_FE_LIBRARIES ${UTOPIA_FE_LIBRARIES} PARENT_SCOPE)
endfunction()

macro(add_utopiafe_app app_name)
    add_library(${app_name} EXCLUDE_FROM_ALL STATIC ${LOCAL_SOURCES})
    target_include_directories(
        ${app_name} BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                  ${CMAKE_CURRENT_BINARY_DIR})
    foreach(app_dependency ${ARGN})
        target_link_libraries(${app_name} PUBLIC ${app_dependency})
    endforeach()
    # foreach(app_dependency ${ARGN}) target_include_directories(${app_name}
    # PUBLIC $<TARGET_PROPERTY:${app_dependency},INTERFACE_INCLUDE_DIRECTORIES>)
    # endforeach()
    set_utopia_compiler_features(${app_name})
    if(NOT TOP_LEVEL)
        set(UTOPIA_FE_APPS
            "${UTOPIA_FE_APPS};${app_name}"
            PARENT_SCOPE)
    endif()
    # install(TARGETS ${app_name} RUNTIME DESTINATION bin LIBRARY DESTINATION
    # lib ARCHIVE DESTINATION lib ) install(FILES ${LOCAL_HEADERS} DESTINATION
    # include)
endmacro(add_utopiafe_app)

###################################

function(append_list_to_string_for_makefile OUTPUT_STRING PRE_LIST POST_LIST
         ITEM_PREFIX LIST_OF_STRINGS)
    # http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro

    list(REMOVE_DUPLICATES LIST_OF_STRINGS)

    set(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

    foreach(LIST_EL ${LIST_OF_STRINGS})

        string(REGEX MATCH "(.*)framework" IS_FRAME_WORK "${LIST_EL}")
        if(IS_FRAME_WORK)
            get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
            string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH
                                 "${LIST_EL}")
            string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME
                                 "${FRAMEWORK_NAME}")

            set(TEMPSTR
                "${TEMPSTR}-F${FRAMEWORK_PATH} -framework ${FRAMEWORK_NAME} ")
        else()
            set(TEMPSTR "${TEMPSTR}${ITEM_PREFIX}${LIST_EL} ")
        endif()
    endforeach(LIST_EL)

    set(TEMPSTR "${TEMPSTR} ${POST_LIST}")
    set(${OUTPUT_STRING}
        ${TEMPSTR}
        PARENT_SCOPE)
endfunction()

function(
    append_list_to_string_for_makefile_with_postfix
    OUTPUT_STRING
    PRE_LIST
    POST_LIST
    ITEM_PREFIX
    ITEM_POSTFIX
    LIST_OF_STRINGS)
    # http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro

    list(REMOVE_DUPLICATES LIST_OF_STRINGS)

    set(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

    foreach(LIST_EL ${LIST_OF_STRINGS})

        string(REGEX MATCH "(.*)framework" IS_FRAME_WORK "${LIST_EL}")
        if(IS_FRAME_WORK)
            get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
            string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH
                                 "${LIST_EL}")
            string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME
                                 "${FRAMEWORK_NAME}")

            set(TEMPSTR
                "${TEMPSTR}-F${FRAMEWORK_PATH} -framework ${FRAMEWORK_NAME} ")
        else()
            set(TEMPSTR "${TEMPSTR}${ITEM_PREFIX}${LIST_EL}${ITEM_POSTFIX} ")
        endif()
    endforeach(LIST_EL)

    set(TEMPSTR "${TEMPSTR} ${POST_LIST}")
    set(${OUTPUT_STRING}
        ${TEMPSTR}
        PARENT_SCOPE)
endfunction()

# ##############################################################################

function(append_list_to_string_for_cmake OUTPUT_STRING PRE_LIST POST_LIST
         LIST_OF_STRINGS)
    # http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro
    list(REMOVE_DUPLICATES LIST_OF_STRINGS)

    set(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

    foreach(LIST_EL ${LIST_OF_STRINGS})
        set(TEMPSTR "${TEMPSTR}\t\"${LIST_EL}\"\n")
    endforeach(LIST_EL)

    set(TEMPSTR "${TEMPSTR} ${POST_LIST}")
    set(${OUTPUT_STRING}
        ${TEMPSTR}
        PARENT_SCOPE)
endfunction()

# ##############################################################################

function(scan_directories in_root_dir in_dirs_to_be_scanned out_includes
         out_headers out_sources)

  # APPEND the results
  set(local_includes ${${out_includes}})
  set(local_headers ${${out_headers}})
  set(local_sources ${${out_sources}})

  find_project_files(${in_root_dir} "${in_dirs_to_be_scanned}" local_headers
                     local_sources)

  foreach(dir ${in_dirs_to_be_scanned})
    list(APPEND local_includes ${in_root_dir}/${dir})
  endforeach(dir)

  set(${out_includes}
      ${local_includes}
      PARENT_SCOPE)
  set(${out_headers}
      ${local_headers}
      PARENT_SCOPE)
  set(${out_sources}
      ${local_sources}
      PARENT_SCOPE)

endfunction()