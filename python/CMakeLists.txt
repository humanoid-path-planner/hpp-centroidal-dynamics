# Define the wrapper library that wraps our library
ADD_LIBRARY(${PY_NAME} SHARED centroidal_dynamics_python.cpp)
TARGET_LINK_LIBRARIES(${PY_NAME} ${PROJECT_NAME} eigenpy::eigenpy)
# don't prepend wrapper library name with lib
SET_TARGET_PROPERTIES(${PY_NAME} PROPERTIES PREFIX "")

IF(APPLE)
  # We need to change the extension for python bindings
  SET_TARGET_PROPERTIES(${PY_NAME} PROPERTIES SUFFIX ".so")
ENDIF(APPLE)

INSTALL(TARGETS ${PY_NAME} DESTINATION ${PYTHON_SITELIB})

ADD_PYTHON_UNIT_TEST("python-centroidal-dynamics" "python/test/binding_tests.py" "python")
