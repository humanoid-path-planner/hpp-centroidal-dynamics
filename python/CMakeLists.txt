# Define the wrapper library that wraps our library
add_library(${PY_NAME} SHARED centroidal_dynamics_python.cpp)
target_link_libraries(${PY_NAME} ${PROJECT_NAME} eigenpy::eigenpy)
# don't prepend wrapper library name with lib
set_target_properties(${PY_NAME} PROPERTIES PREFIX "")

if(APPLE)
  # We need to change the extension for python bindings
  set_target_properties(${PY_NAME} PROPERTIES SUFFIX ".so")
endif(APPLE)

install(TARGETS ${PY_NAME} DESTINATION ${PYTHON_SITELIB})

add_python_unit_test("python-centroidal-dynamics"
                     "python/test/binding_tests.py" "python")
