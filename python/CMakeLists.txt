STRING(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})

ADD_REQUIRED_DEPENDENCY("eigenpy")

# Define the wrapper library that wraps our library
add_library( ${PY_NAME} SHARED centroidal_dynamics_python )
target_link_libraries( ${PY_NAME} ${Boost_LIBRARIES} ${PROJECT_NAME} )
# don't prepend wrapper library name with lib
set_target_properties( ${PY_NAME} PROPERTIES PREFIX "" )

IF(APPLE)
  # We need to change the extension for python bindings
  SET_TARGET_PROPERTIES(${PY_NAME} PROPERTIES SUFFIX ".so")
ENDIF(APPLE)
PKG_CONFIG_USE_DEPENDENCY(${PY_NAME} eigenpy)

INSTALL(
  TARGETS ${PY_NAME} DESTINATION ${PYTHON_SITELIB}
  )

ADD_PYTHON_UNIT_TEST("python-centroidal-dynamics" "python/test/binding_tests.py" "python")
