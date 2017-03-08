cmake_minimum_required( VERSION 2.8 )

include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")

# Define the wrapper library that wraps our library
add_library( centroidal_dynamics SHARED centroidal_dynamics_python.cpp )
target_link_libraries( centroidal_dynamics ${Boost_LIBRARIES} ${PROJECT_NAME} )
# don't prepend wrapper library name with lib
set_target_properties( centroidal_dynamics PROPERTIES PREFIX "" )

#INSTALL(TARGETS centroidal_dynamics DESTINATION lib)
