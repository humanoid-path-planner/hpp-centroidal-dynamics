ADD_UNIT_TEST(static-equilibrium test_static_equilibrium)
PKG_CONFIG_USE_DEPENDENCY(static-equilibrium eigen3)
TARGET_LINK_LIBRARIES(static-equilibrium ${PROJECT_NAME})

ADD_UNIT_TEST(lp-solvers test_LP_solvers)
PKG_CONFIG_USE_DEPENDENCY(lp-solvers eigen3)
TARGET_LINK_LIBRARIES(lp-solvers ${PROJECT_NAME} ${QPOASES_LIBRARY})
