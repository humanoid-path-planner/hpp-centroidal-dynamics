CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# Project properties
SET(PROJECT_NAME hpp-centroidal-dynamics)
SET(PROJECT_DESCRIPTION
  "Utility classes for testing (robust) equilibrium of a system in contact with the environment, and other centroidal dynamics methods."
  )

# Project options
OPTION(BUILD_PYTHON_INTERFACE "Build the python bindings" ON)

# Project configuration
SET(PROJECT_USE_CMAKE_EXPORT TRUE)
SET(CUSTOM_HEADER_DIR "hpp/centroidal-dynamics")
SET(CXX_DISABLE_WERROR TRUE)

# JRL-cmakemodule setup
INCLUDE(cmake/hpp.cmake)
INCLUDE(cmake/boost.cmake)

# Project definition
COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})

IF(BUILD_PYTHON_INTERFACE)
  STRING(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})
  ADD_PROJECT_DEPENDENCY(eigenpy 2.7.12 REQUIRED)
ENDIF(BUILD_PYTHON_INTERFACE)

# Project dependencies
ADD_PROJECT_DEPENDENCY(Eigen3 REQUIRED)

set(CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake/find-external/CDD
  ${PROJECT_SOURCE_DIR}/cmake/find-external/CLP
  ${PROJECT_SOURCE_DIR}/cmake/find-external/qpOASES
  )

find_package(CDD REQUIRED)
find_package(qpOASES REQUIRED)
find_package(CLP)

IF("${CLP_LIBRARY}" STREQUAL "CLP_LIBRARY-NOTFOUND")
  MESSAGE(STATUS "CLP_LIBRARY equal to CLP_LIBRARY-NOTFOUND so I assume CLP was not found ")
ELSE()
  message(STATUS "CLP library found, defining macro CLP_FOUND")
  add_definitions(-DCLP_FOUND)
ENDIF()

# Main Library
SET(${PROJECT_NAME}_HEADERS
  include/${CUSTOM_HEADER_DIR}/local_config.hh
  include/${CUSTOM_HEADER_DIR}/util.hh
  include/${CUSTOM_HEADER_DIR}/logger.hh
  include/${CUSTOM_HEADER_DIR}/solver_LP_abstract.hh
  include/${CUSTOM_HEADER_DIR}/solver_LP_qpoases.hh
  include/${CUSTOM_HEADER_DIR}/solver_LP_clp.hh
  include/${CUSTOM_HEADER_DIR}/centroidal_dynamics.hh
  include/${CUSTOM_HEADER_DIR}/stop-watch.hh
  )

SET(${PROJECT_NAME}_SOURCES
  src/centroidal_dynamics.cpp
  src/solver_LP_abstract.cpp
  src/solver_LP_qpoases.cpp
  src/solver_LP_clp.cpp
  src/util.cpp
  src/logger.cpp
  src/stop-watch.cpp
  )

ADD_LIBRARY(${PROJECT_NAME} SHARED
  ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:include>)
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} SYSTEM PUBLIC
    ${EIGEN3_INCLUDE_DIR} ${CDD_INCLUDE_DIRS} ${qpOASES_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CDD_LIBRARIES} ${qpOASES_LIBRARY})

IF(CLP_FOUND)
  TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} SYSTEM PUBLIC "${CLP_INCLUDE_DIR}")
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CLP_LIBRARIES}
    /usr/lib/libCoinUtils.so)
ENDIF(CLP_FOUND)

INSTALL(TARGETS ${PROJECT_NAME} EXPORT ${TARGETS_EXPORT_NAME} DESTINATION lib)

ADD_SUBDIRECTORY(test)
IF(BUILD_PYTHON_INTERFACE)
  ADD_SUBDIRECTORY(python)
ENDIF(BUILD_PYTHON_INTERFACE)

INSTALL(FILES package.xml DESTINATION share/${PROJECT_NAME})
