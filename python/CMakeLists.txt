cmake_minimum_required( VERSION 2.8 )

include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")

ADD_REQUIRED_DEPENDENCY("eigenpy")


# Define the wrapper library that wraps our library
add_library( centroidal_dynamics SHARED centroidal_dynamics_python.cpp )
target_link_libraries( centroidal_dynamics ${Boost_LIBRARIES} ${PROJECT_NAME} )
# don't prepend wrapper library name with lib
set_target_properties( centroidal_dynamics PROPERTIES PREFIX "" )

IF(APPLE)
	# We need to change the extension for python bindings
	SET_TARGET_PROPERTIES(centroidal_dynamics PROPERTIES SUFFIX ".so")
ENDIF(APPLE)
PKG_CONFIG_USE_DEPENDENCY(centroidal_dynamics eigenpy)

INSTALL(
	FILES ${LIBRARY_OUTPUT_PATH}/centroidal_dynamics.so	DESTINATION ${PYTHON_SITELIB}
)
