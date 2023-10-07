/**
 * flat determinat calculation method test
 * @author Tobias Weber
 * @date jun-2022
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o det det.cpp
 */

#define MATH_USE_FLAT_DET 1

#include "matrix_algos.h"
#include "matrix_conts.h"

using namespace m;
using namespace m_ops;


int main()
{
	using t_real = double;
	using t_vec = vec<t_real, std::vector>;
	using t_mat = mat<t_real, std::vector>;

	t_mat M = m::create<t_mat>({ 1.,-2.,3., -4.,5.,6., 7.,8.,9. });
	t_real d = m::det<t_mat, t_vec>(M);
	std::cout << "M       = " << M << std::endl;
	std::cout << "det     = " << d << std::endl;

	auto [M_inv, ok] = m::inv<t_mat, t_vec>(M);
	t_real d_inv = m::det<t_mat, t_vec>(M_inv);
	std::cout << "M_inv   = " << M_inv << std::endl;
	std::cout << "det_inv = " << d_inv << std::endl;

	return 0;
}
