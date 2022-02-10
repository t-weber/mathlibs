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
	std::cout << "M1 + M2:\n";
	for(std::size_t i=0; i<m3.size1(); ++i)
	{
		for(std::size_t j=0; j<m3.size2(); ++j)
			std::cout << m3(i, j) << " ";
		std::cout << std::endl;
	}

	Matrix<t_real, 2, 2> m4 = m1 * m2;
	std::cout << "M1 * M2:\n";
	for(std::size_t i=0; i<m4.size1(); ++i)
	{
		for(std::size_t j=0; j<m4.size2(); ++j)
			std::cout << m4(i, j) << " ";
		std::cout << std::endl;
	}
}


void test_vector()
{
	std::cout << "\n" << __func__ << std::endl;

	Vector<t_real, 3> v1{}, v2{};
	std::cout << "dim: " << v1.size() << " " << std::endl;

	v1[0] = 1;
	v1[1] = 2;
	v1[2] = 3;

	v2[0] = 7;
	v2[1] = 8;
	v2[2] = 9;

	Vector<t_real, 3> v3 = v1 + v2;
	std::cout << "v1 + v2:\n";
	for(std::size_t i=0; i<v3.size(); ++i)
		std::cout << v3[i] << " ";
	std::cout << std::endl;

	t_real s = v1 * v2;
	std::cout << "v1 * v2 = " << s << std::endl;
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
	std::cout << "M1 + M2:\n";
	for(std::size_t i=0; i<m3.size1(); ++i)
	{
		for(std::size_t j=0; j<m3.size2(); ++j)
			std::cout << m3(i, j) << " ";
		std::cout << std::endl;
	}

	MatrixDyn<t_real> m4 = m1 * m2;
	std::cout << "M1 * M2:\n";
	for(std::size_t i=0; i<m4.size1(); ++i)
	{
		for(std::size_t j=0; j<m4.size2(); ++j)
			std::cout << m4(i, j) << " ";
		std::cout << std::endl;
	}
}


void test_vector_dyn()
{
	std::cout << "\n" << __func__ << std::endl;

	VectorDyn<t_real> v1{3}, v2{3};
	std::cout << "dim: " << v1.size() << " " << std::endl;

	v1[0] = 1;
	v1[1] = 2;
	v1[2] = 3;

	v2[0] = 7;
	v2[1] = 8;
	v2[2] = 9;

	VectorDyn<t_real> v3 = v1 + v2;
	std::cout << "v1 + v2:\n";
	for(std::size_t i=0; i<v3.size(); ++i)
		std::cout << v3[i] << " ";
	std::cout << std::endl;

	t_real s = v1 * v2;
	std::cout << "v1 * v2 = " << s << std::endl;
}


int main()
{
	test_matrix();
	test_matrix_dyn();

	test_vector();
	test_vector_dyn();

	return 0;
}
