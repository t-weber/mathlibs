/**
 * tensor class test
 * @author Tobias Weber
 * @date November 2021
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o tensor tensor.cpp
 */

#include "tensor.h"
#include "tensor_stat.h"

#include <iostream>


using t_real = double;


template<std::size_t ...seq>
void print_seq(const std::index_sequence<seq...>&)
{
	(std::cout << ... << (std::to_string(seq) + " "));
	std::cout << std::endl;
}


void test_tuple()
{
	std::cout << "\n" << __func__ << std::endl;

	auto tup = std::make_tuple(1,2,3,4,5,6,7);
	auto tup2 = remove_from_tuple<1,5>(tup);
	std::cout << "new tuple size: " << std::tuple_size<decltype(tup2)>() << std::endl;
	std::cout << "new tuple: "
		<< std::get<0>(tup2) << " " << std::get<1>(tup2) << " "
		<< std::get<2>(tup2) << " " << std::get<3>(tup2) << " "
		<< std::get<4>(tup2) << std::endl;
}


void test_seq()
{
	std::cout << "\n" << __func__ << std::endl;

	auto seq = seq_cat<std::index_sequence>(
		std::make_index_sequence<5>(), std::make_index_sequence<2>());
	std::cout << "sequence: ";
	print_seq(seq);

	std::cout << "first element: " << seq_first<std::index_sequence>(seq) << std::endl;
	auto seq2 = seq_first<std::index_sequence, 3>(seq);
	std::cout << "first 3 elements: ";
	print_seq(seq2);

	auto seq3 = seq_last<std::index_sequence, 3>(seq);
	std::cout << "last 3 elements: ";
	print_seq(seq3);

	auto seq4 = seq_rm<std::index_sequence, 3>(seq);
	std::cout << "remove element: ";
	print_seq(seq4);
}


void test_tensor()
{
	std::cout << "\n" << __func__ << std::endl;

	m::Tensor<t_real, 2,3> t1{}, t2{};
	std::cout << "dims: " << t1.size<0>() << " " << t1.size<1>() << std::endl;
	std::cout << "total size: " << t1.size() << std::endl;
	std::cout << "order: " << t1.order() << std::endl;

	t1(0,0) = 9;
	t1(0,1) = 8;
	t1(1,1) = 5;

	t2(0, 1) = -1;

	t1 += 2.*t2;

	std::cout << "elements: ";
	for(std::size_t i=0; i<t1.size(); ++i)
		std::cout << t1[i] << " ";
	std::cout << std::endl;


	auto t3 = t1 * t2;
	std::cout << "tensor product order: " << t3.order() << std::endl;


	m::Tensor<t_real, 3,3> t4{};
	t4(0,0) = 1; t4(0,1) = 2; t4(0,2) = 3;
	t4(1,0) = 4; t4(1,1) = 5; t4(1,2) = 6;
	t4(2,0) = 7; t4(2,1) = 8; t4(2,2) = 9;
	auto c = t4.contract<0,1>();
	std::cout << "contraction: order: " << c.order() << ", size: " << c.size() << ", value: " << c[0] << std::endl;

	std::cout << "static element access: " << t4.operator()<1,1>() << std::endl;
	std::cout << "dynamic element access: " << t4(1,1) << std::endl;
}


void test_tensor_dyn()
{
	std::cout << "\n" << __func__ << std::endl;

	m::TensorDyn<t_real> t1{{2,3}}, t2{{2,3}};
	std::cout << "dims: " << t1.size(0) << " " << t1.size(1) << std::endl;
	std::cout << "total size: " << t1.size() << std::endl;
	std::cout << "order: " << t1.order() << std::endl;

	t1({0,0}) = 9;
	t1({0,1}) = 8;
	t1({1,1}) = 5;

	t2({0, 1}) = -1;

	t1 += 2.*t2;

	std::cout << "elements: ";
	for(std::size_t i=0; i<t1.size(); ++i)
		std::cout << t1[i] << " ";
	std::cout << std::endl;


	auto t3 = t1 * t2;
	std::cout << "tensor product order: " << t3.order() << std::endl;


	m::TensorDyn<t_real> t4({3,3});
	t4({0,0}) = 1; t4({0,1}) = 2; t4({0,2}) = 3;
	t4({1,0}) = 4; t4({1,1}) = 5; t4({1,2}) = 6;
	t4({2,0}) = 7; t4({2,1}) = 8; t4({2,2}) = 9;

	auto c = t4.contract(0, 1);

	std::cout << "contraction: order: " << c.order() << ", size: " << c.size() << ", value: " << c[0] << std::endl;
	std::cout << "element access: " << t4({1,1}) << std::endl;
}


int main()
{
	test_tuple();
	test_seq();
	test_tensor();
	test_tensor_dyn();

	return 0;
}
