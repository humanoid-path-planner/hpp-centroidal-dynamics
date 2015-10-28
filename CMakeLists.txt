cmake_minimum_required(VERSION 2.6)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/cpack.cmake)


SET(PROJECT_NAME robust-equilibrium-lib)
SET(PROJECT_DESCRIPTION
   "Utility classes for testing (robust) static equilibrium of a system in contact with the environment"
 )
SET(PROJECT_URL "")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${PROJECT_SOURCE_DIR}/cmake2/")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")

set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

set(EXTERNAL_LIBRARY_DIR "${PROJECT_SOURCE_DIR}/external/lib")

if ( MSVC )  
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX")
endif()

SETUP_PROJECT()

string (REPLACE "-Werror" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )

# Declare Headers
SET(${PROJECT_NAME}_HEADERS
    include/robust-equilibrium-lib/config.hh
    include/robust-equilibrium-lib/util.hh
    include/robust-equilibrium-lib/logger.hh
    include/robust-equilibrium-lib/solver_LP_abstract.hh
    include/robust-equilibrium-lib/solver_LP_qpoases.hh
    include/robust-equilibrium-lib/solver_LP_clp.hh
    include/robust-equilibrium-lib/static_equilibrium.hh
  )

find_package(Eigen3)
if (EIGEN3_INCLUDE_DIR)
        message(STATUS "Found Eigen3")
else()
        message(STATUS "Eigen3 not found looking for Eigen2")
        find_package(Eigen2 REQUIRED)
endif()

find_package(CDD REQUIRED)
find_package(CLP REQUIRED)

add_subdirectory (src)
add_subdirectory (test)

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()

