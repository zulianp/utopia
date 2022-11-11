# 2020-07-17 AF: No need to make a target for rapidxml, even less need to install it. We use it only in implementation files, so it's a PRIVATE dependency

# include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#
# list(APPEND LOCAL_HEADERS
#     rapidxml.hpp
#     rapidxml_iterators.hpp
#     rapidxml_print.hpp
#     rapidxml_utils.hpp)
#
#
# install(FILES ${LOCAL_HEADERS} DESTINATION include)
