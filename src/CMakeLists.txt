cmake_minimum_required(VERSION 2.6)

if(CLP_FOUND)
  include_directories("${CLP_INCLUDE_DIR}")
endif()


SET(LIBRARY_NAME ${PROJECT_NAME})

SET(${LIBRARY_NAME}_SOURCES
  centroidal_dynamics.cpp
  solver_LP_abstract.cpp
  solver_LP_qpoases.cpp
  solver_LP_clp.cpp
  util.cpp
  logger.cpp
  stop-watch.cpp
  )

ADD_LIBRARY(${LIBRARY_NAME} SHARED ${${LIBRARY_NAME}_SOURCES})

TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${CDD_LIBRARIES} ${qpOASES_LIBRARY})
PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} eigen3)

if(CLP_FOUND)
  TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${CLP_LIBRARIES}
    /usr/lib/libCoinUtils.so)
endif()

INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION lib)
