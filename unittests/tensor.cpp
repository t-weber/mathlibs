/**
 * tensor class test
 * @author Tobias Weber
 * @date November 2021
 * @license see 'LICENSE' file
 *
 * References:
 *  * https://www.boost.org/doc/libs/1_76_0/libs/test/doc/html/index.html
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o tensor tensor.cpp
 */

#define BOOST_TEST_MODULE test_tensor

#include <iostream>

#include <boost/test/included/unit_test.hpp>
#include <boost/type_index.hpp>
namespace test = boost::unit_test;
namespace ty = boost::typeindex;

#include "tensor.h"
#include "tensor_stat.h"


BOOST_AUTO_TEST_CASE_TEMPLATE(test_tensor, t_real, decltype(std::tuple<float, double, long double>{}))
{
	t_real eps = std::pow(std::numeric_limits<t_real>::epsilon(), 0.5);
	std::cout << "\nRunning \"test_tensor\" with "
		<< ty::type_id_with_cvr<t_real>().pretty_name()
		<< " type and epsilon = " << eps << "." << std::endl;

	m::Tensor<t_real, 2, 3> t1{}, t2{};
	std::cout << "dims: " << t1.template size<0>() << " " << t1.template size<1>() << std::endl;
	std::cout << "total size: " << t1.size() << std::endl;
	std::cout << "order: " << t1.order() << std::endl;

	t1(0, 0) = 9;
	t1(0, 1) = 8;
	t1(1, 1) = 5;

	t2(0, 1) = -1;

	t1 += t_real{2} * t2;

	std::cout << "elements: ";
	for(std::size_t i = 0; i < t1.size(); ++i)
		std::cout << t1[i] << " ";
	std::cout << std::endl;


	auto t3 = t1 * t2;
	std::cout << "tensor product order: " << t3.order() << std::endl;


	m::Tensor<t_real, 3, 3> t4{};
	t4(0,0) = 1; t4(0,1) = 2; t4(0,2) = 3;
	t4(1,0) = 4; t4(1,1) = 5; t4(1,2) = 6;
	t4(2,0) = 7; t4(2,1) = 8; t4(2,2) = 9;
	auto c = t4.template contract<0, 1>();
	std::cout << "contraction: order: " << c.order()
		<< ", size: " << c.size() << ", value: " << c[0]
		<< std::endl;
	BOOST_TEST((c.order() == 0));
	BOOST_TEST((c.size() == 1));
	BOOST_TEST((m::equals<t_real>(c[0], t_real{15}, eps)));


	std::cout << "static element access: " << t4.template operator()<1, 1>() << std::endl;
	std::cout << "dynamic element access: " << t4(1, 1) << std::endl;
	BOOST_TEST((m::equals<t_real>(t4.template operator()<1, 1>(), t_real{5}, eps)));
	BOOST_TEST((m::equals<t_real>(t4(1, 1), t_real{5}, eps)));
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_tensor_dyn, t_real, decltype(std::tuple<float, double, long double>{}))
{
	t_real eps = std::pow(std::numeric_limits<t_real>::epsilon(), 0.5);
	std::cout << "\nRunning \"test_tensor_dyn\" with "
		<< ty::type_id_with_cvr<t_real>().pretty_name()
		<< " type and epsilon = " << eps << "." << std::endl;

	m::TensorDyn<t_real> t1{{2, 3}}, t2{{2, 3}};
	std::cout << "dims: " << t1.size(0) << " " << t1.size(1) << std::endl;
	std::cout << "total size: " << t1.size() << std::endl;
	std::cout << "order: " << t1.order() << std::endl;

	t1({0, 0}) = 9;
	t1({0, 1}) = 8;
	t1({1, 1}) = 5;

	t2({0, 1}) = -1;

	t1 += t_real{2} * t2;

	std::cout << "elements: ";
	for(std::size_t i = 0; i < t1.size(); ++i)
		std::cout << t1[i] << " ";
	std::cout << std::endl;


	auto t3 = t1 * t2;
	std::cout << "tensor product order: " << t3.order() << std::endl;


	m::TensorDyn<t_real> t4({3, 3});
	t4({0, 0}) = 1; t4({0, 1}) = 2; t4({0, 2}) = 3;
	t4({1, 0}) = 4; t4({1, 1}) = 5; t4({1, 2}) = 6;
	t4({2, 0}) = 7; t4({2, 1}) = 8; t4({2, 2}) = 9;

	auto c = t4.contract(0, 1);
	std::cout << "contraction: order: " << c.order() << ", size: " << c.size() << ", value: " << c[0] << std::endl;
	BOOST_TEST((c.order() == 0));
	BOOST_TEST((c.size() == 1));
	BOOST_TEST((m::equals<t_real>(c[0], t_real{15}, eps)));


	std::cout << "element access: " << t4({1,1}) << std::endl;
	BOOST_TEST((m::equals<t_real>(t4({1, 1}), t_real{5}, eps)));
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_tensor2, t_real, decltype(std::tuple<float, double, long double>{}))
{
	t_real eps = std::pow(std::numeric_limits<t_real>::epsilon(), 0.5);
	std::cout << "\nRunning \"test_tensor2\" with "
		<< ty::type_id_with_cvr<t_real>().pretty_name()
		<< " type and epsilon = " << eps << "." << std::endl;

	m::Tensor<t_real, 2, 2, 3> t1{}, t2{};
	std::cout << "dims: " << t1.template size<0>() << " " << t1.template size<1>()
		<< " " << t1.template size<2>() << std::endl;
	std::cout << "total size: " << t1.size() << std::endl;
	std::cout << "order: " << t1.order() << std::endl;

	t1(1, 0, 0) = 9;
	t1(1, 0, 1) = 8;
	t1(1, 1, 1) = 5;

	t2(1, 0, 1) = -1;

	t1 += t_real{2} * t2;

	std::cout << "elements: ";
	for(std::size_t i = 0; i < t1.size(); ++i)
		std::cout << t1[i] << " ";
	std::cout << std::endl;


	auto t3 = t1 * t2;
	std::cout << "tensor product order: " << t3.order() << std::endl;


	m::Tensor<t_real, 2, 3, 3> t4{};
	t4(1, 0, 0) = 1; t4(1, 0, 1) = 2; t4(1, 0, 2) = 3;
	t4(1, 1, 0) = 4; t4(1, 1, 1) = 5; t4(1, 1, 2) = 6;
	t4(1, 2, 0) = 7; t4(1, 2, 1) = 8; t4(1, 2, 2) = 9;
	auto c = t4.template contract<1, 2>();
	std::cout << "contraction: order: " << c.order()
		<< ", size: " << c.size() << ", values: " << c[0] << " " << c[1]
		<< std::endl;
	BOOST_TEST((c.order() == 1));
	BOOST_TEST((c.size() == 2));
	BOOST_TEST((m::equals<t_real>(c[0], t_real{0}, eps)));
	BOOST_TEST((m::equals<t_real>(c[1], t_real{15}, eps)));


	std::cout << "static element access: " << t4.template operator()<1, 1, 1>() << std::endl;
	std::cout << "dynamic element access: " << t4(1, 1, 1) << std::endl;
	BOOST_TEST((m::equals<t_real>(t4.template operator()<1, 1, 1>(), t_real{5}, eps)));
	BOOST_TEST((m::equals<t_real>(t4(1, 1, 1), t_real{5}, eps)));
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_tensor_dyn2, t_real, decltype(std::tuple<float, double, long double>{}))
{
	t_real eps = std::pow(std::numeric_limits<t_real>::epsilon(), 0.5);
	std::cout << "\nRunning \"test_tensor_dyn2\" with "
		<< ty::type_id_with_cvr<t_real>().pretty_name()
		<< " type and epsilon = " << eps << "." << std::endl;

	m::TensorDyn<t_real> t1{{2, 2, 3}}, t2{{2, 2, 3}};
	std::cout << "dims: " << t1.size(0) << " " << t1.size(1) << " " << t1.size(2) << std::endl;
	std::cout << "total size: " << t1.size() << std::endl;
	std::cout << "order: " << t1.order() << std::endl;

	t1({1, 0, 0}) = 9;
	t1({1, 0, 1}) = 8;
	t1({1, 1, 1}) = 5;

	t2({1, 0, 1}) = -1;

	t1 += t_real{2} * t2;

	std::cout << "elements: ";
	for(std::size_t i = 0; i < t1.size(); ++i)
		std::cout << t1[i] << " ";
	std::cout << std::endl;


	auto t3 = t1 * t2;
	std::cout << "tensor product order: " << t3.order() << std::endl;


	m::TensorDyn<t_real> t4({2, 3, 3});
	t4({1, 0, 0}) = 1; t4({1, 0, 1}) = 2; t4({1, 0, 2}) = 3;
	t4({1, 1, 0}) = 4; t4({1, 1, 1}) = 5; t4({1, 1, 2}) = 6;
	t4({1, 2, 0}) = 7; t4({1, 2, 1}) = 8; t4({1, 2, 2}) = 9;

	auto c = t4.contract(1, 2);

	std::cout << "contraction: order: " << c.order()
		<< ", size: " << c.size() << ", values: " << c[0] << " " << c[1] << std::endl;
	BOOST_TEST((c.order() == 1));
	BOOST_TEST((c.size() == 2));
	BOOST_TEST((m::equals<t_real>(c[0], t_real{0}, eps)));
	BOOST_TEST((m::equals<t_real>(c[1], t_real{15}, eps)));


	std::cout << "element access: " << t4({1, 1, 1}) << std::endl;
	BOOST_TEST((m::equals<t_real>(t4({1, 1, 1}), t_real{5}, eps)));
}
