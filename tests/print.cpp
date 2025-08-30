/**
 * printing test
 * @author Tobias Weber
 * @date 30-aug-2025
 * @license see 'LICENSE' file
 *
 * g++ -std=c++23 -I../libs -Wall -Wextra -Weffc++ -o print print.cpp
 */

#include "matrix_algos.h"
#include "matrix_conts.h"
using namespace m_ops;

#if __cplusplus >= 202302L
	#include <print>
#endif


int main()
{
	using t_real = double;
	using t_vec = m::vec<t_real, std::vector>;
	using t_mat = m::mat<t_real, std::vector>;

	const t_real pi = m::pi<t_real>;
	t_vec v = m::create<t_vec>({ pi, 2.*pi, 3.*pi, 4.*pi });
	t_mat M = m::B_matrix<t_mat>(4., 5., 6., 90./180.*pi, 90./180.*pi, 60./180.*pi);

	std::cout.precision(4);
	std::cout << "vector: " << v << std::endl;
	std::cout << "matrix: " << M << std::endl;

#if __cplusplus >= 202302L
	std::println("vector: {:.4f}", v);
	std::println("matrix: {:.4f}", M);
#endif
	return 0;
}
