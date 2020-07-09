function(format_setup)
  if(NOT CLANG_FORMAT_EXE)
    set(CLANG_FORMAT_EXE clang-format)
  endif()

  if(NOT EXISTS ${CLANG_FORMAT_EXE})
    find_program(format_executable_tmp NAMES ${CLANG_FORMAT_EXE} clang-format-mp-7.0 clang-format-9  clang-format-mp-10 HINTS /opt/local/bin/ /usr/local/bin/ /usr/bin/ )
    if(format_executable_tmp)
      set(CLANG_FORMAT_EXE ${format_executable_tmp})
      unset(format_executable_tmp)
    else()
      message(WARNING "ClangFormat: ${CLANG_FORMAT_EXE} not found!")
      return()
    endif()
  endif()

  foreach(format_source ${ARGV})
    get_filename_component(format_source ${format_source} ABSOLUTE)
    list(APPEND format_sources ${format_source})
  endforeach()

  add_custom_target(utopia_format
    COMMAND
      ${CLANG_FORMAT_EXE}
      -style=file -fallback-style=none
      -i
      ${format_sources}
    WORKING_DIRECTORY
      ${CMAKE_SOURCE_DIR}
    COMMENT
      "Formating with ${CLANG_FORMAT_EXE} ..."
  )

  if(TARGET format)
    add_dependencies(format utopia_format)
  else()
    add_custom_target(format DEPENDS utopia_format)
  endif()
endfunction()

function(target_format_setup target)
  get_target_property(target_sources ${target} SOURCES)
  format_setup(${target_sources})
endfunction()
