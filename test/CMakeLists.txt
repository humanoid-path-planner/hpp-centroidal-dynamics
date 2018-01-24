cmake_minimum_required(VERSION 2.6)


INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${CDD_INCLUDE_DIR}")


# Make Boost.Test generates the main function in test cases.
#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

# ADD_TESTCASE(NAME)
# ------------------------
#
# Define a test named `NAME'.
#
# This macro will create a binary from `NAME.cc', link it against
# Boost and add it to the test suite.
#
MACRO(ADD_TESTCASE NAME GENERATED)
  IF (${GENERATED} STREQUAL TRUE)
    ADD_EXECUTABLE(${NAME} ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cpp)
  ELSE()
    ADD_EXECUTABLE(${NAME} ${NAME}.cpp)
  ENDIF()
  ADD_TEST(${NAME} ${RUNTIME_OUTPUT_DIRECTORY}/${NAME})
 
  # Link against Boost and project library.
  TARGET_LINK_LIBRARIES(${NAME}
    ${Boost_LIBRARIES}
    ${PROJECT_NAME}
    )

ENDMACRO(ADD_TESTCASE)


if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )

ADD_TESTCASE(test_static_equilibrium FALSE)
ADD_TESTCASE(test_LP_solvers FALSE)

