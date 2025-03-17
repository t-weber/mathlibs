/**
 * sequence test
 * @author Tobias Weber
 * @date November 2021
 * @license see 'LICENSE' file
 *
 * References:
 *  * https://www.boost.org/doc/libs/1_76_0/libs/test/doc/html/index.html
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o seq seq.cpp
 */

#define BOOST_TEST_MODULE test_tensor

#include <iostream>

#include <boost/test/included/unit_test.hpp>

#include "variadic_algos.h"


template<std::size_t ...seq>
void print_seq(const std::index_sequence<seq...>&)
{
	(std::cout << ... << (std::to_string(seq) + " "));
	std::cout << std::endl;
}


BOOST_AUTO_TEST_CASE(test_tuple)
{
	std::cout << "\n" << __func__ << std::endl;

	auto tup = std::make_tuple(1, 2, 3, 4, 5, 6, 7);
	auto tup2 = remove_from_tuple<1, 5>(tup);
	std::cout << "new tuple size: " << std::tuple_size<decltype(tup2)>() << std::endl;
	std::cout << "new tuple: "
		<< std::get<0>(tup2) << " " << std::get<1>(tup2) << " "
		<< std::get<2>(tup2) << " " << std::get<3>(tup2) << " "
		<< std::get<4>(tup2) << std::endl;

	BOOST_TEST((std::tuple_size<decltype(tup2)>() == 5));
}


BOOST_AUTO_TEST_CASE(test_seq)
{
	std::cout << "\n" << __func__ << std::endl;

	auto seq = seq_cat<std::index_sequence>(
		std::make_index_sequence<5>(), std::make_index_sequence<2>());
	std::cout << "sequence: ";
	print_seq(seq);

	std::cout << "first element: " << seq_first_elem<std::index_sequence>(seq) << std::endl;
	auto seq2 = seq_first<std::index_sequence, 3>(seq);
	std::cout << "first 3 elements: ";
	print_seq(seq2);

	auto seq3 = seq_last<std::index_sequence, 3>(seq);
	std::cout << "last 3 elements: ";
	print_seq(seq3);

	auto seq4 = seq_rm<std::index_sequence, 3>(seq);
	std::cout << "remove element: ";
	print_seq(seq4);

	BOOST_TEST((seq_first_elem<std::index_sequence>(seq) == 0));
}
