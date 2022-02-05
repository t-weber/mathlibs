/**
 * matrix class test
 * @author Tobias Weber
 * @date Februray 2022
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o matrix matrix.cpp
 */

#include "tensor.h"
#include "tensor_dyn.h"

#include <iostream>


using t_real = double;


void test_matrix()
{
	std::cout << "\n" << __func__ << std::endl;

	Matrix<t_real, 2, 2> m1{}, m2{};
	std::cout << "dims: " << m1.size1() << " " << m1.size2() << std::endl;

	m1(0,0) = 1; m1(0,1) = 2;
	m1(1,0) = 3; m1(1,1) = 4;

	m2(0,0) = 9; m2(0,1) = 8;
	m2(1,0) = 7; m2(1,1) = 6;

	Matrix<t_real, 2, 2> m3 = m1 + m2;

	std::cout << "elements:\n";
	for(std::size_t i=0; i<m3.size1(); ++i)
	{
		for(std::size_t j=0; j<m3.size2(); ++j)
			std::cout << m3(i, j) << " ";
		std::cout << std::endl;
	}
}


void test_matrix_dyn()
{
	std::cout << "\n" << __func__ << std::endl;

	MatrixDyn<t_real> m1{2, 2}, m2{2, 2};
	std::cout << "dims: " << m1.size1() << " " << m1.size2() << std::endl;

	m1(0,0) = 1; m1(0,1) = 2;
	m1(1,0) = 3; m1(1,1) = 4;

	m2(0,0) = 9; m2(0,1) = 8;
	m2(1,0) = 7; m2(1,1) = 6;

	MatrixDyn<t_real> m3 = m1 + m2;

	std::cout << "elements:\n";
	for(std::size_t i=0; i<m3.size1(); ++i)
	{
		for(std::size_t j=0; j<m3.size2(); ++j)
			std::cout << m3(i, j) << " ";
		std::cout << std::endl;
	}
}


int main()
{
	test_matrix();
	test_matrix_dyn();

	return 0;
}
