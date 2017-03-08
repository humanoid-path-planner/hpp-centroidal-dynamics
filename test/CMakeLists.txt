cmake_minimum_required(VERSION 2.6)

include_directories("${SRC_DIR}")
include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${CDD_INCLUDE_DIR}")
# no need to include directories for qpOASES as it is automatically done through pkgconfig
if(CLP_FOUND)
  include_directories("${CLP_INCLUDE_DIR}")
endif()

PROJECT(centroidal-dynamics-lib)

if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )

add_executable(test_static_equilibrium test_static_equilibrium.cpp)
add_executable(test_LP_solvers test_LP_solvers.cpp)

TARGET_LINK_LIBRARIES(test_LP_solvers centroidal-dynamics-lib)
TARGET_LINK_LIBRARIES(test_static_equilibrium centroidal-dynamics-lib)
#~ TARGET_LINK_LIBRARIES(polytopetest polytope ${SRC_DIR}/../external/cddlib-094b/lib-src/libcdd.a)
