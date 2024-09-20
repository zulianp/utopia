# -*- mode: cmake -*-
#
# UTOPIA Find Module
#
# Usage:
#    Control the search through UTOPIA_DIR.
#
#    Following variables are set:
#    UTOPIA_FOUND            (BOOL)               Flag indicating if UTOPIA was found
#    UTOPIA_INCLUDES         (LIST of PATH)       Path to the UTOPIA include file
#    UTOPIA_DEFS             (LIST)               List with all defs for the compiler
#    UTOPIA_LIBRARIES        (LIST)               List of all required UTOPIA libraries
#    Utopia::utopia          TARGET               Exported target to be used in target_link_libararies, no need to fiddle with UTOPIA_INCLUDES/UTOPIA_DEFS/UTOPIA_LIBRARIES
#
# #############################################################################

set(UTOPIA_SEARCH_PATHS "/usr/local/;/usr/;${Utopia_DIR};${UTOPIA_DIR}
        ${UTOPIA_INCLUDES}")

if(UTOPIA_ENABLE_ENV_READ)
  set(UTOPIA_SEARCH_PATHS "${UTOPIA_SEARCH_PATHS};$ENV{UTOPIA_DIR};$ENV{Utopia_DIR}")
endif()

find_package(Utopia CONFIG
    HINTS ${UTOPIA_SEARCH_PATHS}
    PATH_SUFFIXES lib/cmake
)

