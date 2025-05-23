# ##############################################################################
# ##############################################################################
# ##############################################################################
# Utitlity functions

# ##############################################################################

# Place for utopia_add_library

# ##############################################################################

function(create_absolute_paths root paths absolutePaths)
  foreach(path ${${paths}})
    set(temp; ${temp} ${${root}}/${path})
  endforeach()
  set(${absolutePaths}
      ${temp}
      PARENT_SCOPE)
endfunction()

# ##############################################################################

function(find_project_files rootPath dirPaths headers sources)
  set(verbose TRUE)

  set(theaders ${${headers}})
  set(tsources ${${sources}})

  foreach(INCLUDE_PATH ${dirPaths})
    # INCLUDE_DIRECTORIES(${rootPath}/${INCLUDE_PATH})

    file(GLOB TEMP_HPPSRC CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.cpp")
    file(GLOB TEMP_SRC CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.c")
    file(GLOB TEMP_HPPDR CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.hpp")
    file(GLOB TEMP_HDR CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.h")

    source_group(${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR};
                                       ${TEMP_HPPSRC}; ${TEMP_SRC}; ${TEMP_UI})

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

function(append_list_to_string_for_makefile OUTPUT_STRING PRE_LIST POST_LIST
         ITEM_PREFIX LIST_OF_STRINGS)
  # http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro
  list(REMOVE_DUPLICATES LIST_OF_STRINGS)

  set(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

  foreach(LIST_EL ${LIST_OF_STRINGS})

    string(REGEX MATCH "(.*)framework" IS_FRAME_WORK "${LIST_EL}")
    if(IS_FRAME_WORK)
      get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
      string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH "${LIST_EL}")
      string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME "${FRAMEWORK_NAME}")
      set(TEMPSTR "${TEMPSTR}-F${FRAMEWORK_PATH} -framework ${FRAMEWORK_NAME} ")
    else()
      set(TEMPSTR "${TEMPSTR}${ITEM_PREFIX}${LIST_EL} ")
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

# ##############################################################################

function (getListOfVarsStartingWith _prefix _varResult)
    get_cmake_property(_vars VARIABLES)
    string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    set (${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

# ##############################################################################

function(string_starts_with str search result)
  string(FIND "${str}" "${search}" out)

  if("${out}" EQUAL 0)
    set(${result} TRUE)
    set(${result} TRUE PARENT_SCOPE)
  else()
    set(${result} FALSE PARENT_SCOPE)
  endif()
endfunction()