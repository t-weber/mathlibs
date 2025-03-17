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

	add_executable(utst_seq unittests/seq.cpp)
	target_link_libraries(utst_seq mathlibs)
	add_test(utst_seq utst_seq)

	add_executable(utst_tensor unittests/tensor.cpp)
	target_link_libraries(utst_tensor mathlibs)
	add_test(utst_tensor utst_tensor)

	add_executable(utst_matrix unittests/matrix.cpp)
	target_link_libraries(utst_matrix mathlibs)
	add_test(utst_matrix utst_matrix)

	add_executable(utst_rotation unittests/rotation.cpp)
	target_link_libraries(utst_rotation mathlibs)
	add_test(utst_rotation utst_rotation)
endif()
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# test programs
# -----------------------------------------------------------------------------
if(BUILD_TESTS)
	add_executable(tst_hash tests/hash.cpp)
	target_link_libraries(tst_hash mathlibs)

	add_executable(tst_matrix tests/matrix.cpp)
	target_link_libraries(tst_matrix mathlibs)

	add_executable(tst_vec tests/vec.cpp)
	target_link_libraries(tst_vec mathlibs)

	add_executable(tst_quat tests/quat.cpp)
	target_link_libraries(tst_quat mathlibs)

	add_executable(tst_det tests/det.cpp)
	target_link_libraries(tst_det mathlibs)

	add_executable(tst_intersect tests/intersect.cpp)
	target_link_libraries(tst_intersect mathlibs)

	add_executable(tst_phys tests/phys.cpp)
	target_link_libraries(tst_phys mathlibs)

	add_executable(tst_recursion tests/recursion.cpp)
	target_link_libraries(tst_recursion mathlibs)

	#add_executable(tst_triag tests/triag.cpp)
	#target_link_libraries(tst_triag mathlibs)
endif()
# -----------------------------------------------------------------------------
