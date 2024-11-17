/**
 * physics algos tests
 * @author Tobias Weber
 * @date 3 November 2024
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o phys phys.cpp
 */

#include "tensor.h"
#include "tensor_stat.h"
#include "matrix_algos.h"
#include "matrix_conts.h"
#include "phys_algos.h"

#include <complex>
#include <iostream>



/**
 * test of topology functions using eigenvector matrix
 */
void test_topo()
{
	using namespace m_ops;
	std::cout << "\n" << __func__ << std::endl;

	using t_real = double;
	using t_cplx = std::complex<t_real>;
	using t_mat = m::Matrix<t_cplx, 2, 2>;
	using t_mat3 = m::Matrix<t_cplx, 3, 3>;
	using t_vec = m::Vector<t_cplx, 3>;
	using t_vec_real = m::Vector<t_real, 3>;

	auto get_state = [](const t_vec_real& /*Q*/) -> t_mat
	{
		// TODO
		t_mat state{};
		state(0,0) = 1; state(0,1) = 2;
		state(1,0) = 3; state(1,1) = 4;
		return state;
	};

	t_vec_real Q = m::create<t_vec_real>({ 0, 0, 0 });
	std::vector<t_vec> conns =
		m::berry_connection<t_mat, t_vec, t_vec_real>(get_state, Q, 0.01);
	std::cout << "conn = [ ";
	for(std::size_t i = 0; i < conns[0].size(); ++i)
		std::cout << conns[0][i] << " ";
	std::cout << "]" << std::endl;

	std::vector<t_cplx> curvs_2d =
		m::berry_curvature_2d<t_mat, t_vec, t_vec_real>(get_state, Q, 0.01);
	std::cout << "curv = [ ";
	for(const t_cplx& curv : curvs_2d)
		std::cout << curv << " ";
	std::cout << "]" << std::endl;

	std::vector<t_mat3> curvs =
		m::berry_curvature<t_mat, t_mat3, t_vec, t_vec_real>(get_state, Q, 0.01);
	std::cout << "curs =\n";
	for(const t_mat3& curv : curvs)
		std::cout << curv << std::endl;

	std::vector<t_cplx> nums_2d =
		m::chern_numbers_2d<t_mat, t_vec, t_vec_real>(get_state, Q, 0.5, 0.01, 0.01);
	std::cout << "chern (via boundary) = [ ";
		for(const t_cplx& curv : nums_2d)
			std::cout << curv << " ";
	std::cout << "]" << std::endl;

	std::vector<t_cplx> nums_2d_ar =
		m::chern_numbers_2d_area<t_mat, t_vec, t_vec_real>(get_state, Q, 0.5, 0.01, 0.01);
	std::cout << "chern (via area) = [ ";
		for(const t_cplx& curv : nums_2d_ar)
			std::cout << curv << " ";
	std::cout << "]" << std::endl;

	std::cout << std::endl;
}



/**
 * test of topology functions using orthonormal eigenvectors
 */
void test_topo2()
{
	using namespace m_ops;
	std::cout << "\n" << __func__ << std::endl;

	using t_real = double;
	using t_cplx = std::complex<t_real>;
	using t_mat = m::Matrix<t_cplx, 3, 3>;
	using t_vec = m::Vector<t_cplx, 3>;
	using t_vec_real = m::Vector<t_real, 3>;

	auto get_state = [](const t_vec_real& /*Q*/) -> std::vector<t_vec>
	{
		// TODO
		std::vector<t_vec> states{};

		states.emplace_back(m::create<t_vec>({ 1, 3 }));
		states.emplace_back(m::create<t_vec>({ 2, 4 }));

		return states;
	};

	t_vec_real Q = m::create<t_vec_real>({ 0, 0, 0 });
	std::vector<t_vec> conns =
		m::berry_connection<t_vec, t_vec_real>(get_state, Q, 0.01);
	std::cout << "conn = [ ";
	for(std::size_t i = 0; i < conns[0].size(); ++i)
		std::cout << conns[0][i] << " ";
	std::cout << "]" << std::endl;

	std::vector<t_cplx> curvs_2d =
		m::berry_curvature_2d<t_vec, t_vec_real>(get_state, Q, 0.01);
	std::cout << "curv = [ ";
	for(const t_cplx& curv : curvs_2d)
		std::cout << curv << " ";
	std::cout << "]" << std::endl;

	std::vector<t_mat> curvs =
		m::berry_curvature<t_mat, t_vec, t_vec_real>(get_state, Q, 0.01);
	std::cout << "curs =\n";
	for(const t_mat& curv : curvs)
		std::cout << curv << std::endl;

	std::cout << std::endl;
}



int main()
{
	try
	{
		test_topo();
		test_topo2();
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return -1;
	}

	return 0;
}
