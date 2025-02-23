#
# mathlibs
# @author Tobias Weber (orcid: 0000-0002-7230-1932)
# @date 22-feb-2025
# @license see 'LICENSE' file
#

cmake_minimum_required(VERSION 3.5)
#cmake_policy(SET CMP0167 NEW)


# -----------------------------------------------------------------------------
# project options
# -----------------------------------------------------------------------------
project(mathlibs)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)


option(BUILD_UNIT_TESTS "build unit tests"    TRUE)
option(BUILD_TESTS      "build test programs" TRUE)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# dependencies
# -----------------------------------------------------------------------------
find_package(Boost REQUIRED)
add_definitions(${Boost_CXX_FLAGS})


include(GNUInstallDirs)

include_directories(
	${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/libs
	${CMAKE_CURRENT_BINARY_DIR}
	${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}/..
)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# mathlibs
# -----------------------------------------------------------------------------
add_library(mathlibs INTERFACE)

target_sources(mathlibs
	PUBLIC
		${PROJECT_SOURCE_DIR}/libs/tensor.h
		${PROJECT_SOURCE_DIR}/libs/tensor_stat.h
		${PROJECT_SOURCE_DIR}/libs/variadic_algos.h

		${PROJECT_SOURCE_DIR}/libs/matrix_concepts.h
		${PROJECT_SOURCE_DIR}/libs/matrix_conts.h
		${PROJECT_SOURCE_DIR}/libs/matrix_algos.h

		${PROJECT_SOURCE_DIR}/libs/graph_conts.h
		${PROJECT_SOURCE_DIR}/libs/graph_algos.h

		${PROJECT_SOURCE_DIR}/libs/geo_conts.h
		${PROJECT_SOURCE_DIR}/libs/geo_algos.h

		${PROJECT_SOURCE_DIR}/libs/phys_algos.h
		${PROJECT_SOURCE_DIR}/libs/poly_algos.h
		${PROJECT_SOURCE_DIR}/libs/str_algos.h
		${PROJECT_SOURCE_DIR}/libs/tree_algos.h

		${PROJECT_SOURCE_DIR}/libs/helpers.h

	PRIVATE
)


set_target_properties(mathlibs PROPERTIES PUBLIC_HEADER
	"${PROJECT_SOURCE_DIR}/libs/tensor.h;\
${PROJECT_SOURCE_DIR}/libs/tensor_stat.h;\
${PROJECT_SOURCE_DIR}/libs/variadic_algos.h;\
${PROJECT_SOURCE_DIR}/libs/matrix_concepts.h;\
${PROJECT_SOURCE_DIR}/libs/matrix_conts.h;\
${PROJECT_SOURCE_DIR}/libs/matrix_algos.h;\
${PROJECT_SOURCE_DIR}/libs/graph_conts.h;\
${PROJECT_SOURCE_DIR}/libs/graph_algos.h;\
${PROJECT_SOURCE_DIR}/libs/geo_conts.h;\
${PROJECT_SOURCE_DIR}/libs/geo_algos.h;\
${PROJECT_SOURCE_DIR}/libs/phys_algos.h;\
${PROJECT_SOURCE_DIR}/libs/poly_algos.h;\
${PROJECT_SOURCE_DIR}/libs/str_algos.h;\
${PROJECT_SOURCE_DIR}/libs/tree_algos.h;\
${PROJECT_SOURCE_DIR}/libs/helpers.h"
)


install(TARGETS mathlibs
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mathlibs
)


# cmake file
install(FILES ${PROJECT_SOURCE_DIR}/FindMathlibs.cmake
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake
)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# unit tests
# -----------------------------------------------------------------------------
if(BUILD_UNIT_TESTS)
	enable_testing()

	add_executable(matrix1 unittests/matrix1.cpp)
	target_link_libraries(matrix1 mathlibs)
	add_test(matrix1 matrix1)

	add_executable(rotation unittests/rotation.cpp)
	target_link_libraries(rotation mathlibs)
	add_test(rotation rotation)
endif()
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# test programs
# -----------------------------------------------------------------------------
if(BUILD_TESTS)
	add_executable(hash tests/hash.cpp)
	target_link_libraries(hash mathlibs)

	add_executable(tensor tests/tensor.cpp)
	target_link_libraries(tensor mathlibs)

	add_executable(matrix tests/matrix.cpp)
	target_link_libraries(matrix mathlibs)

	add_executable(vec tests/vec.cpp)
	target_link_libraries(vec mathlibs)

	add_executable(quat tests/quat.cpp)
	target_link_libraries(quat mathlibs)

	add_executable(det tests/det.cpp)
	target_link_libraries(det mathlibs)

	add_executable(intersect tests/intersect.cpp)
	target_link_libraries(intersect mathlibs)

	add_executable(phys tests/phys.cpp)
	target_link_libraries(phys mathlibs)

	add_executable(recursion tests/recursion.cpp)
	target_link_libraries(recursion mathlibs)

	#add_executable(triag tests/triag.cpp)
	#target_link_libraries(triag mathlibs)
endif()
# -----------------------------------------------------------------------------
