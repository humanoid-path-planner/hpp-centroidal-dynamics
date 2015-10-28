# - Try to find libcdd
# Once done this will define
#  CDD_FOUND - System has CDD
#  CDD_INCLUDE_DIRS - The CDD include directories
#  CDD_LIBRARIES - The libraries needed to use CDD
#  CDD_DEFINITIONS - Compiler switches required for using CDD

# /usr/include/coin, /usr/lib/libClp.so

find_path(CLP_INCLUDE_DIR coin/ClpSimplex.hpp
          HINTS ${CLP_INCLUDEDIR}
          PATH_SUFFIXES CLP )

find_library(CLP_LIBRARY NAMES libclp
             HINTS ${CLP_LIBDIR} ${CLP_LIBRARY_DIRS} )

set(CLP_LIBRARIES ${CLP_LIBRARY} )
set(CLP_INCLUDE_DIRS ${CLP_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CDD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CLP  DEFAULT_MSG
                                  CLP_LIBRARY CLP_INCLUDE_DIR)

mark_as_advanced(CLP_INCLUDE_DIR CLP_LIBRARY )
