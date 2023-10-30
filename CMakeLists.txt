cmake_minimum_required(VERSION 3.10)

# Project properties
set(PROJECT_NAME hpp-centroidal-dynamics)
set(PROJECT_DESCRIPTION
    "Utility classes for testing (robust) equilibrium of a system in contact with the environment, and other centroidal dynamics methods."
)

# Project options
option(BUILD_PYTHON_INTERFACE "Build the python bindings" ON)

# Project configuration
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CUSTOM_HEADER_DIR "hpp/centroidal-dynamics")
set(CXX_DISABLE_WERROR TRUE)

# JRL-cmakemodule setup
include(cmake/hpp.cmake)
include(cmake/boost.cmake)

# Project definition
compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})

if(BUILD_PYTHON_INTERFACE)
  string(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})
  add_project_dependency(eigenpy 2.7.12 REQUIRED)
endif(BUILD_PYTHON_INTERFACE)

# Project dependencies
add_project_dependency(Eigen3 REQUIRED)

set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/find-external/CDD
    ${PROJECT_SOURCE_DIR}/cmake/find-external/CLP
    ${PROJECT_SOURCE_DIR}/cmake/find-external/qpOASES)

find_package(CDD REQUIRED)
find_package(qpOASES REQUIRED)
find_package(CLP)

if("${CLP_LIBRARY}" STREQUAL "CLP_LIBRARY-NOTFOUND")
  message(
    STATUS
      "CLP_LIBRARY equal to CLP_LIBRARY-NOTFOUND so I assume CLP was not found "
  )
else()
  message(STATUS "CLP library found, defining macro CLP_FOUND")
  add_definitions(-DCLP_FOUND)
endif()

# Main Library
set(${PROJECT_NAME}_HEADERS
    include/${CUSTOM_HEADER_DIR}/local_config.hh
    include/${CUSTOM_HEADER_DIR}/util.hh
    include/${CUSTOM_HEADER_DIR}/logger.hh
    include/${CUSTOM_HEADER_DIR}/solver_LP_abstract.hh
    include/${CUSTOM_HEADER_DIR}/solver_LP_qpoases.hh
    include/${CUSTOM_HEADER_DIR}/solver_LP_clp.hh
    include/${CUSTOM_HEADER_DIR}/centroidal_dynamics.hh
    include/${CUSTOM_HEADER_DIR}/stop-watch.hh)

set(${PROJECT_NAME}_SOURCES
    src/centroidal_dynamics.cpp
    src/solver_LP_abstract.cpp
    src/solver_LP_qpoases.cpp
    src/solver_LP_clp.cpp
    src/util.cpp
    src/logger.cpp
    src/stop-watch.cpp)

add_library(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SOURCES}
                                   ${${PROJECT_NAME}_HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:include>)
target_include_directories(
  ${PROJECT_NAME} SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR} ${CDD_INCLUDE_DIRS}
                                ${qpOASES_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${CDD_LIBRARIES} ${qpOASES_LIBRARY})

if(CLP_FOUND)
  target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC "${CLP_INCLUDE_DIR}")
  target_link_libraries(${PROJECT_NAME} ${CLP_LIBRARIES}
                        /usr/lib/libCoinUtils.so)
endif(CLP_FOUND)

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${TARGETS_EXPORT_NAME}
  DESTINATION lib)

add_subdirectory(test)
if(BUILD_PYTHON_INTERFACE)
  add_subdirectory(python)
endif(BUILD_PYTHON_INTERFACE)

install(FILES package.xml DESTINATION share/${PROJECT_NAME})
