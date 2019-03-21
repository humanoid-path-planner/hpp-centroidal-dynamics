cmake_minimum_required(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/test.cmake)
INCLUDE(cmake/python.cmake)
INCLUDE(cmake/hpp.cmake)

SET(PROJECT_NAME hpp-centroidal-dynamics)
SET(PROJECT_DESCRIPTION
    "Utility classes for testing (robust) equilibrium of a system in contact with the environment, and other centroidal dynamics methods."
    )

SET(CUSTOM_HEADER_DIR "hpp/centroidal-dynamics")

SETUP_HPP_PROJECT()

# Inhibit all warning messages.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# remove flag that makes all warnings into errors
string (REPLACE "-Werror" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )

OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" ON)
IF(BUILD_PYTHON_INTERFACE)
    # search for python
    FINDPYTHON(2.7 EXACT REQUIRED)
    find_package( PythonLibs 2.7 REQUIRED )
    include_directories( ${PYTHON_INCLUDE_DIRS} )

    ADD_REQUIRED_DEPENDENCY("eigenpy")

    find_package( Boost COMPONENTS python REQUIRED )
    include_directories( ${Boost_INCLUDE_DIR} )

ENDIF(BUILD_PYTHON_INTERFACE)

ADD_REQUIRED_DEPENDENCY("eigen3")
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/find-external/CDD
    ${PROJECT_SOURCE_DIR}/cmake/find-external/CLP
    ${PROJECT_SOURCE_DIR}/cmake/find-external/qpOASES
    )

find_package(CDD REQUIRED)
find_package(qpOASES REQUIRED)
find_package(CLP)

IF("${CLP_LIBRARY}" STREQUAL "CLP_LIBRARY-NOTFOUND")
    message(STATUS "CLP_LIBRARY equal to CLP_LIBRARY-NOTFOUND so I assume CLP was not found ")
else()
    message(STATUS "CLP library found, defining macro CLP_FOUND")
    add_definitions(-DCLP_FOUND)
endif()

INCLUDE_DIRECTORIES(${qpOASES_INCLUDE_DIRS})


ADD_SUBDIRECTORY (include/${CUSTOM_HEADER_DIR})
add_subdirectory (src)
add_subdirectory (test)
IF(BUILD_PYTHON_INTERFACE)
    add_subdirectory (python)
ENDIF(BUILD_PYTHON_INTERFACE)

SETUP_PROJECT_FINALIZE()
