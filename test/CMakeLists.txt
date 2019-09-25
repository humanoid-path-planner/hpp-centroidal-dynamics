cmake_minimum_required(VERSION 2.6)


# Make Boost.Test generates the main function in test cases.
#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

ADD_UNIT_TEST(static-equilibrium test_static_equilibrium)
PKG_CONFIG_USE_DEPENDENCY(static-equilibrium eigen3)
TARGET_LINK_LIBRARIES(static-equilibrium ${PROJECT_NAME})

ADD_UNIT_TEST(lp-solvers test_LP_solvers)
PKG_CONFIG_USE_DEPENDENCY(lp-solvers eigen3)
TARGET_LINK_LIBRARIES(lp-solvers ${PROJECT_NAME} ${QPOASES_LIBRARY})

IF(NOT ${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    ADD_CUSTOM_TARGET(link_target ALL COMMAND ${CMAKE_COMMAND} -E create_symlink
        "${PROJECT_SOURCE_DIR}/test_data"
        "${PROJECT_BINARY_DIR}/test_data")
ENDIF(NOT ${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
