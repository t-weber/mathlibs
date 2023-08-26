/**
 * vector tests
 * @author Tobias Weber
 * @date aug-2023
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o vec vec.cpp
 */

#include <vector>
#include <iostream>

#include "math_algos.h"
#include "math_conts.h"

using namespace m;
using namespace m_ops;


int main()
{
	using t_real = double;
	using t_vec = vec<t_real, std::vector>;
	//using t_mat = mat<t_real, std::vector>;

	std::vector<t_vec> vecs
	{{
		m::create<t_vec>({ 1., 2., 3. }),
		m::create<t_vec>({ 5., 6., 7. }),
		m::create<t_vec>({ -9., 8., -7. }),
		m::create<t_vec>({ -9., 1., 4. }),
	}};

	auto [min, max] = minmax_comp(vecs);
	std::cout << "min: " << min << "\nmax: " << max << std::endl;

	t_vec mid = avg<t_vec>(vecs);
	auto [min_dist, max_dist] = minmax_dist(vecs, mid);
	std::cout << "mid: " << mid << "\nmin_dist: " << min_dist << ", max_dist: " << max_dist << std::endl;

	return 0;
}
