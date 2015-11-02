cmake_minimum_required(VERSION 2.6)

include_directories("${SRC_DIR}")
include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${CDD_INCLUDE_DIR}")
include_directories("${qpOASES_INCLUDE_DIR}")

if( DEFINED CLP_FOUND)
  include_directories("${CLP_INCLUDE_DIR}")
endif()


SET(LIBRARY_NAME ${PROJECT_NAME})

SET(${LIBRARY_NAME}_SOURCES
    ${INCLUDE_DIR}/robust-equilibrium-lib/config.hh
    ${INCLUDE_DIR}/robust-equilibrium-lib/solver_LP_abstract.hh
    ${INCLUDE_DIR}/robust-equilibrium-lib/solver_LP_qpoases.hh
    ${INCLUDE_DIR}/robust-equilibrium-lib/solver_LP_clp.hh
    ${INCLUDE_DIR}/robust-equilibrium-lib/static_equilibrium.hh
    static_equilibrium.cpp
    solver_LP_abstract.cpp
    solver_LP_qpoases.cpp
    solver_LP_clp.cpp
    util.cpp
    logger.cpp
    stop-watch.cpp
  )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
ADD_LIBRARY(${LIBRARY_NAME} SHARED ${${LIBRARY_NAME}_SOURCES})

if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )

TARGET_LINK_LIBRARIES(robust-equilibrium-lib ${CDD_LIBRARY}
                                             ${qpOASES_LIBRARY})
if( DEFINED CLP_FOUND)
  TARGET_LINK_LIBRARIES(robust-equilibrium-lib ${CLP_LIBRARY}
                                               /usr/lib/libCoinUtils.so)
endif()

SET_TARGET_PROPERTIES(robust-equilibrium-lib PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET_TARGET_PROPERTIES(robust-equilibrium-lib PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION lib)
