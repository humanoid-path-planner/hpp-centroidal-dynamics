add_unit_test(static-equilibrium test_static_equilibrium.cpp)
target_link_libraries(static-equilibrium ${PROJECT_NAME})

add_unit_test(lp-solvers test_LP_solvers.cpp)
target_link_libraries(lp-solvers ${PROJECT_NAME})

if(NOT ${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
  add_custom_target(
    link_target ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${PROJECT_SOURCE_DIR}/test_data" "${PROJECT_BINARY_DIR}/test_data")
endif(NOT ${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
