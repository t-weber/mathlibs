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
#include "phys_algos.h"

#include <complex>
#include <iostream>


void test_topo()
{
	std::cout << "\n" << __func__ << std::endl;

	using t_real = double;
	using t_cplx = std::complex<t_real>;
	using t_mat = m::Matrix<t_cplx, 2, 2>;
	using t_mat_real = m::Matrix<t_real, 2, 2>;
	using t_vec = m::Vector<t_cplx, 3>;
	using t_vec_real = m::Vector<t_real, 3>;

	auto get_state = [](const t_vec_real& Q) -> t_mat
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

	std::vector<t_cplx> curvs =
		m::berry_curvature<t_mat, t_vec, t_vec_real>(get_state, Q, 0.01);
	std::cout << "curv = [ ";
	for(const t_cplx& curv : curvs)
		std::cout << curv << " ";
	std::cout << "]" << std::endl;

	std::cout << std::endl;
}


int main()
{
	try
	{
		test_topo();
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return -1;
	}

	return 0;
}
