/**
 * rotation matrix unit tests
 * @author Tobias Weber
 * @date march-2022
 * @license see 'LICENSE' file
 *
 * References:
 *  * https://www.boost.org/doc/libs/1_76_0/libs/test/doc/html/index.html
 *
 * g++ -std=c++20 -I .. -o rotation rotation.cpp
 */

#define BOOST_TEST_MODULE test_rotation

#include <iostream>
#include <tuple>

#include <boost/test/included/unit_test.hpp>
#include <boost/type_index.hpp>
namespace test = boost::unit_test;
namespace ty = boost::typeindex;

#include "libs/math_algos.h"
#include "libs/math_conts.h"


BOOST_AUTO_TEST_CASE_TEMPLATE(test_rotation,
	t_scalar, decltype(std::tuple<float, double, long double>{}))
{
	using namespace m_ops;

	t_scalar eps = std::pow(std::numeric_limits<t_scalar>::epsilon(), 0.5);

	std::cout << "Testing with " << ty::type_id_with_cvr<t_scalar>().pretty_name()
		<< " type and epsilon = " << eps << "." << std::endl;

	using t_vector = m::vec<t_scalar, std::vector>;
	using t_matrix = m::mat<t_scalar, std::vector>;


	t_vector vecFrom = m::create<t_vector>({1, 0, 1});
	t_vector vecTo = m::create<t_vector>({1, 0, 0});

	t_matrix mat_1a = m::rotation<t_matrix, t_vector, t_scalar>(vecFrom, vecTo, eps);
	t_matrix mat_1b = m::rotation_nd<t_matrix, t_vector, t_scalar>(vecFrom, vecTo);

	BOOST_TEST((m::equals<t_matrix>(mat_1a, mat_1b, eps)));
	std::cout << mat_1a << std::endl;
	std::cout << mat_1b << std::endl;


	t_vector vecFrom2 = m::create<t_vector>({1, 0, 1, 0, 1});
	t_vector vecTo2 = m::create<t_vector>({1, 0, 0, 0, 0});
	t_matrix mat_2 = m::rotation_nd<t_matrix, t_vector, t_scalar>(vecFrom2, vecTo2);

	t_vector vecTo2b = (mat_2 * vecFrom2);
	vecTo2b *= m::norm<t_vector>(vecTo2) / m::norm<t_vector>(vecFrom2);
	std::cout << "vec: " << vecTo2b << std::endl;
	BOOST_TEST((m::equals<t_vector>(vecTo2, vecTo2b, eps)));
}
