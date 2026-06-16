cmake_minimum_required(VERSION 3.22)

# Project properties
set(PROJECT_NAME hpp-centroidal-dynamics)
set(PROJECT_DESCRIPTION
    "Utility classes for testing (robust) equilibrium of a system in contact with the environment, and other centroidal dynamics methods."
)

# Project configuration
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CUSTOM_HEADER_DIR "hpp/centroidal-dynamics")
set(CXX_DISABLE_WERROR TRUE)

# Check if the submodule cmake have been initialized
set(JRL_CMAKE_MODULES "${CMAKE_CURRENT_LIST_DIR}/cmake")
if(EXISTS "${JRL_CMAKE_MODULES}/base.cmake")
    message(STATUS "JRL cmakemodules found in 'cmake/' git submodule")
else()
    find_package(jrl-cmakemodules QUIET CONFIG)
    if(jrl-cmakemodules_FOUND)
        get_property(
            JRL_CMAKE_MODULES
            TARGET jrl-cmakemodules::jrl-cmakemodules
            PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        )
        message(
            STATUS
            "JRL cmakemodules found on system at ${JRL_CMAKE_MODULES}"
        )
    else()
        message(STATUS "JRL cmakemodules not found. Let's fetch it.")
        include(FetchContent)
        FetchContent_Declare(
            "jrl-cmakemodules"
            GIT_REPOSITORY "https://github.com/jrl-umi3218/jrl-cmakemodules.git"
        )
        FetchContent_MakeAvailable("jrl-cmakemodules")
        FetchContent_GetProperties(
            "jrl-cmakemodules"
            SOURCE_DIR JRL_CMAKE_MODULES
        )
    endif()
endif()

include("${JRL_CMAKE_MODULES}/hpp.cmake")
include("${JRL_CMAKE_MODULES}/boost.cmake")

# Project definition
compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})

if(BUILD_PYTHON_INTERFACE)
    string(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})
    add_project_dependency(eigenpy 2.7.12 REQUIRED)
endif(BUILD_PYTHON_INTERFACE)

# Project dependencies
add_project_dependency(Eigen3 REQUIRED)
add_project_dependency(Boost REQUIRED)

set(CMAKE_MODULE_PATH
    ${JRL_CMAKE_MODULES}/find-external/CDD
    ${JRL_CMAKE_MODULES}/find-external/CLP
    ${JRL_CMAKE_MODULES}/find-external/CoinUtils
    ${JRL_CMAKE_MODULES}/find-external/qpOASES
)

add_project_dependency(CDD REQUIRED FIND_EXTERNAL CDD)
add_project_dependency(CLP QUIET FIND_EXTERNAL CLP)
add_project_dependency(qpOASES QUIET FIND_EXTERNAL qpOASES)

if(CLP_FOUND)
    add_project_dependency(CoinUtils REQUIRED FIND_EXTERNAL CoinUtils)
    message(STATUS "Building with coin-or clp")
else()
    message(STATUS "Building without coin-or clp")
endif()

if(qpOASES_FOUND)
    message(STATUS "Building with qpOASES")
else()
    message(STATUS "Building without qpOASES")
endif()

if(NOT CLP_FOUND AND NOT qpOASES_FOUND)
    message(FATAL_ERROR "Need CLP or qpOASES")
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
    include/${CUSTOM_HEADER_DIR}/stop-watch.hh
)

set(${PROJECT_NAME}_SOURCES
    src/centroidal_dynamics.cpp
    src/solver_LP_abstract.cpp
    src/solver_LP_qpoases.cpp
    src/solver_LP_clp.cpp
    src/util.cpp
    src/logger.cpp
    src/stop-watch.cpp
)

add_library(
    ${PROJECT_NAME}
    SHARED
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_HEADERS}
)
target_include_directories(
    ${PROJECT_NAME}
    PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC Eigen3::Eigen CDD::CDD Boost::boost
)

if(CLP_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC CoinUtils::CoinUtils CLP::CLP)
    target_compile_definitions(${PROJECT_NAME} PUBLIC CLP_FOUND)
endif()

if(qpOASES_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC qpOASES::qpOASES)
    target_compile_definitions(${PROJECT_NAME} PUBLIC qpOASES_FOUND)
endif()

install(TARGETS ${PROJECT_NAME} EXPORT ${TARGETS_EXPORT_NAME} DESTINATION lib)

add_subdirectory(test)
if(BUILD_PYTHON_INTERFACE)
    add_subdirectory(python)
endif(BUILD_PYTHON_INTERFACE)
