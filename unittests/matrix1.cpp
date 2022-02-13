/**
 * matrix unit tests
 * @author Tobias Weber
 * @date february-2022
 * @license see 'LICENSE' file
 *
 * References:
 *  * https://www.boost.org/doc/libs/1_76_0/libs/test/doc/html/index.html
 *
 * g++ -std=c++20 -I .. -o matrix1 matrix1.cpp
 */

#define BOOST_TEST_MODULE test_matrix1

#include <iostream>
#include <tuple>

#include <boost/test/included/unit_test.hpp>
#include <boost/type_index.hpp>
namespace test = boost::unit_test;
namespace ty = boost::typeindex;

#include "libs/tensor.h"
#include "libs/tensor_stat.h"
#include "libs/math_algos.h"


BOOST_AUTO_TEST_CASE_TEMPLATE(test_matrix1, t_scalar, decltype(std::tuple<float, double, long double>{}))
{
	t_scalar eps = std::pow(std::numeric_limits<t_scalar>::epsilon(), 0.5);

	std::cout << "Testing with " << ty::type_id_with_cvr<t_scalar>().pretty_name()
		<< " type and epsilon = " << eps << "." << std::endl;

	using t_matrix = m::MatrixDyn<t_scalar>;
	using t_matrix_static_33 = m::Matrix<t_scalar, 3, 3>;

	t_matrix mat1(3, 3);
	t_matrix_static_33 mat1s;

	// LinearAlgebra.inv([1 2 2; 5 9 5; 3 7 7])
	mat1(0,0)=mat1s(0,0)=1.; mat1(0,1)=mat1s(0,1)=2.; mat1(0,2)=mat1s(0,2)=2.;
	mat1(1,0)=mat1s(1,0)=5.; mat1(1,1)=mat1s(1,1)=9.; mat1(1,2)=mat1s(1,2)=5.;
	mat1(2,0)=mat1s(2,0)=3.; mat1(2,1)=mat1s(2,1)=7.; mat1(2,2)=mat1s(2,2)=7.;

	t_matrix mat1_inv = mat1.inverse();
	t_matrix_static_33 mat1s_inv = mat1s.inverse();

	BOOST_TEST((mat1_inv.size1() == mat1s_inv.size1()));
	BOOST_TEST((mat1_inv.size2() == mat1s_inv.size2()));

	for(std::size_t i=0; i<mat1_inv.size1(); ++i)
	{
		for(std::size_t j=0; j<mat1_inv.size2(); ++j)
		{
			BOOST_TEST((m::equals<t_scalar>(mat1_inv(i, j), mat1s_inv(i, j), eps)));
		}
	}

	BOOST_TEST((m::equals<t_scalar>(mat1_inv(0,0), 7., eps)));
	BOOST_TEST((m::equals<t_scalar>(mat1_inv(0,1), 0., eps)));
	BOOST_TEST((m::equals<t_scalar>(mat1_inv(0,2), -2., eps)));
	BOOST_TEST((m::equals<t_scalar>(mat1_inv(1,0), -5., eps)));
	BOOST_TEST((m::equals<t_scalar>(mat1_inv(1,1), 0.25, eps)));
	BOOST_TEST((m::equals<t_scalar>(mat1_inv(1,2), 1.25, eps)));
	BOOST_TEST((m::equals<t_scalar>(mat1_inv(2,0), 2., eps)));
	BOOST_TEST((m::equals<t_scalar>(mat1_inv(2,1), -0.25, eps)));
	BOOST_TEST((m::equals<t_scalar>(mat1_inv(2,2), -0.25, eps)));
}
