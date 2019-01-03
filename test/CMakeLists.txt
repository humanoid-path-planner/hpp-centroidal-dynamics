cmake_minimum_required(VERSION 2.6)


# Make Boost.Test generates the main function in test cases.
#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

ADD_UNIT_TEST(static-equilibrium test_static_equilibrium)
PKG_CONFIG_USE_DEPENDENCY(static-equilibrium eigen3)
TARGET_LINK_LIBRARIES(static-equilibrium ${PROJECT_NAME})

ADD_UNIT_TEST(lp-solvers test_LP_solvers)
PKG_CONFIG_USE_DEPENDENCY(lp-solvers qpOASES)
PKG_CONFIG_USE_DEPENDENCY(lp-solvers eigen3)
TARGET_LINK_LIBRARIES(lp-solvers ${PROJECT_NAME})
