ADD_UNIT_TEST(static-equilibrium test_static_equilibrium.cpp)
TARGET_LINK_LIBRARIES(static-equilibrium ${PROJECT_NAME})

ADD_UNIT_TEST(lp-solvers test_LP_solvers.cpp)
TARGET_LINK_LIBRARIES(lp-solvers ${PROJECT_NAME})

IF(NOT ${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    ADD_CUSTOM_TARGET(link_target ALL COMMAND ${CMAKE_COMMAND} -E create_symlink
        "${PROJECT_SOURCE_DIR}/test_data"
        "${PROJECT_BINARY_DIR}/test_data")
ENDIF(NOT ${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
