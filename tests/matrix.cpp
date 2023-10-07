/**
 * (tensor-class derived) matrix class tests
 * @author Tobias Weber
 * @date Februray 2022
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o matrix matrix.cpp
 */

#include "tensor.h"
#include "tensor_stat.h"
#include "matrix_algos.h"

#include <iostream>


using t_real = double;


void test_matrix()
{
	std::cout << "\n" << __func__ << std::endl;

	using t_mat22 = m::Matrix<t_real, 2, 2>;
	using t_mat11 = m::Matrix<t_real, 1, 1>;
	using t_vec2 = m::Vector<t_real, 2>;

	t_mat22 m1{}, m2{};
	std::cout << "dims: " << m1.size1() << " " << m1.size2() << std::endl;

	m1(0,0) = 1; m1(0,1) = 2;
	m1(1,0) = 3; m1(1,1) = 4;

	m2(0,0) = 9; m2(0,1) = 8;
	m2(1,0) = 7; m2(1,1) = 6;

	t_mat22 m3 = m1 + m2;
	std::cout << "M1 + M2:\n";
	for(std::size_t i=0; i<m3.size1(); ++i)
	{
		for(std::size_t j=0; j<m3.size2(); ++j)
			std::cout << m3(i, j) << " ";
		std::cout << std::endl;
	}

	t_mat22 m4 = m1 * m2;
	std::cout << "M1 * M2:\n";
	for(std::size_t i=0; i<m4.size1(); ++i)
	{
		for(std::size_t j=0; j<m4.size2(); ++j)
			std::cout << m4(i, j) << " ";
		std::cout << std::endl;
	}

	t_mat11 sub = m1.submatrix<0,1>();
	std::cout << "submatrix:\n";
	for(std::size_t i=0; i<sub.size1(); ++i)
	{
		for(std::size_t j=0; j<sub.size2(); ++j)
			std::cout << sub(i, j) << " ";
		std::cout << std::endl;
	}

	std::cout << "determinant: " << m1.determinant() << std::endl;
	std::cout << "determinant (general function): " << m::det<t_mat22, t_vec2>(m1) << std::endl;

	t_mat22 invmat = m1.inverse();
	std::cout << "inverse:\n";
	for(std::size_t i=0; i<invmat.size1(); ++i)
	{
		for(std::size_t j=0; j<invmat.size2(); ++j)
			std::cout << invmat(i, j) << " ";
		std::cout << std::endl;
	}

	/*auto [invmat2, invmat2_ok] = m::inv<t_mat22, t_vec2>(m1);
	std::cout << "inverse (general function), ok = " << invmat2_ok << ":\n";
	for(std::size_t i=0; i<invmat2.size1(); ++i)
	{
		for(std::size_t j=0; j<invmat2.size2(); ++j)
			std::cout << invmat2(i, j) << " ";
		std::cout << std::endl;
	}*/

	t_mat22 ident = invmat * m1;
	std::cout << "identity:\n";
	for(std::size_t i=0; i<ident.size1(); ++i)
	{
		for(std::size_t j=0; j<ident.size2(); ++j)
			std::cout << ident(i, j) << " ";
		std::cout << std::endl;
	}

	t_mat22 transp = m1.transpose();
	std::cout << "transposed:\n";
	for(std::size_t i=0; i<transp.size1(); ++i)
	{
		for(std::size_t j=0; j<transp.size2(); ++j)
			std::cout << transp(i, j) << " ";
		std::cout << std::endl;
	}

	t_mat22 transp2 = m::trans(m1);
	std::cout << "transposed (general function):\n";
	for(std::size_t i=0; i<transp2.size1(); ++i)
	{
		for(std::size_t j=0; j<transp2.size2(); ++j)
			std::cout << transp2(i, j) << " ";
		std::cout << std::endl;
	}
}


void test_vector()
{
	std::cout << "\n" << __func__ << std::endl;

	m::Vector<t_real, 3> v1{}, v2{};
	std::cout << "dim: " << v1.size() << " " << std::endl;

	v1[0] = 1;
	v1[1] = 2;
	v1[2] = 3;

	v2[0] = 7;
	v2[1] = 8;
	v2[2] = 9;

	m::Vector<t_real, 3> v3 = v1 + v2;
	std::cout << "v1 + v2:\n";
	for(std::size_t i=0; i<v3.size(); ++i)
		std::cout << v3[i] << " ";
	std::cout << std::endl;

	t_real s = v1 * v2;
	std::cout << "v1 * v2 = " << s << std::endl;
}


void test_row_col_vector()
{
	std::cout << "\n" << __func__ << std::endl;

	m::ColVector<t_real, 3> v1{};
	m::RowVector<t_real, 3> v2{};
	std::cout << "dim: " << v1.size() << " " << std::endl;

	v1[0] = 1;
	v1[1] = 2;
	v1[2] = 3;

	v2[0] = 7;
	v2[1] = 8;
	v2[2] = 9;

	m::ColVector<t_real, 3> v3 = v1 + v1;
	std::cout << "v1 + v1:\n";
	for(std::size_t i=0; i<v3.size(); ++i)
		std::cout << v3[i] << " ";
	std::cout << std::endl;

	t_real s = v1 * v2;
	std::cout << "v1 * v2 = " << s << std::endl;

	v2 = v1.transpose();
	for(std::size_t i=0; i<v2.size(); ++i)
		std::cout << v2[i] << " ";
	std::cout << std::endl;
}


void test_matrix_dyn()
{
	std::cout << "\n" << __func__ << std::endl;

	using t_mat = m::MatrixDyn<t_real>;
	//using t_vec = m::VectorDyn<t_real>;

	t_mat m1{2, 2}, m2{2, 2};
	std::cout << "dims: " << m1.size1() << " " << m1.size2() << std::endl;

	m1(0,0) = 1; m1(0,1) = 2;
	m1(1,0) = 3; m1(1,1) = 4;

	m2(0,0) = 9; m2(0,1) = 8;
	m2(1,0) = 7; m2(1,1) = 6;

	t_mat m3 = m1 + m2;
	std::cout << "M1 + M2:\n";
	for(std::size_t i=0; i<m3.size1(); ++i)
	{
		for(std::size_t j=0; j<m3.size2(); ++j)
			std::cout << m3(i, j) << " ";
		std::cout << std::endl;
	}

	t_mat m4 = m1 * m2;
	std::cout << "M1 * M2:\n";
	for(std::size_t i=0; i<m4.size1(); ++i)
	{
		for(std::size_t j=0; j<m4.size2(); ++j)
			std::cout << m4(i, j) << " ";
		std::cout << std::endl;
	}

	t_mat sub = m1.submatrix(0,1);
	std::cout << "submatrix:\n";
	for(std::size_t i=0; i<sub.size1(); ++i)
	{
		for(std::size_t j=0; j<sub.size2(); ++j)
			std::cout << sub(i, j) << " ";
		std::cout << std::endl;
	}

	std::cout << "determinant: " << m1.determinant() << std::endl;

	t_mat invmat = m1.inverse();
	std::cout << "inverse:\n";
	for(std::size_t i=0; i<invmat.size1(); ++i)
	{
		for(std::size_t j=0; j<invmat.size2(); ++j)
			std::cout << invmat(i, j) << " ";
		std::cout << std::endl;
	}

	/*auto [invmat2, invmat2_ok] = m::inv<t_mat, t_vec>(m1);
	std::cout << "inverse (direct call), ok = " << invmat2_ok << ":\n";
	for(std::size_t i=0; i<invmat2.size1(); ++i)
	{
		for(std::size_t j=0; j<invmat2.size2(); ++j)
			std::cout << invmat2(i, j) << " ";
		std::cout << std::endl;
	}*/

	t_mat ident = invmat * m1;
	std::cout << "identity:\n";
	for(std::size_t i=0; i<ident.size1(); ++i)
	{
		for(std::size_t j=0; j<ident.size2(); ++j)
			std::cout << ident(i, j) << " ";
		std::cout << std::endl;
	}
}


void test_vector_dyn()
{
	std::cout << "\n" << __func__ << std::endl;

	m::VectorDyn<t_real> v1{3}, v2{3};
	std::cout << "dim: " << v1.size() << " " << std::endl;

	v1[0] = 1;
	v1[1] = 2;
	v1[2] = 3;

	v2[0] = 7;
	v2[1] = 8;
	v2[2] = 9;

	m::VectorDyn<t_real> v3 = v1 + v2;
	std::cout << "v1 + v2:\n";
	for(std::size_t i=0; i<v3.size(); ++i)
		std::cout << v3[i] << " ";
	std::cout << std::endl;

	t_real s = v1 * v2;
	std::cout << "v1 * v2 = " << s << std::endl;
}


void test_row_col_vector_dyn()
{
	std::cout << "\n" << __func__ << std::endl;

	m::ColVectorDyn<t_real> v1{3};
	m::RowVectorDyn<t_real> v2{3};
	std::cout << "dim: " << v1.size() << " " << std::endl;

	v1[0] = 1;
	v1[1] = 2;
	v1[2] = 3;

	v2[0] = 7;
	v2[1] = 8;
	v2[2] = 9;

	m::ColVectorDyn<t_real> v3 = v1 + v1;
	std::cout << "v1 + v1:\n";
	for(std::size_t i=0; i<v3.size(); ++i)
		std::cout << v3[i] << " ";
	std::cout << std::endl;

	t_real s = v1 * v2;
	std::cout << "v1 * v2 = " << s << std::endl;

	v2 = v1.transpose();
	for(std::size_t i=0; i<v2.size(); ++i)
		std::cout << v2[i] << " ";
	std::cout << std::endl;
}


int main()
{
	try
	{
		test_matrix();
		test_matrix_dyn();

		test_vector();
		test_vector_dyn();

		test_row_col_vector();
		test_row_col_vector_dyn();
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return -1;
	}

	return 0;
}
