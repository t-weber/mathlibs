/**
 * quaternion test
 * @author Tobias Weber
 * @date may-2022
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o intersect intersect.cpp
 */

#include "matrix_algos.h"
#include "matrix_conts.h"

using namespace m;
using namespace m_ops;


template<class t_scalar, class t_vec>
void intersect_line_plane()
{
	t_scalar eps = 1e-6;

	t_vec lineOrg = m::create<t_vec>({ 1, 0, 0, });
	t_vec lineDir = m::create<t_vec>({ 0, 0, 5, });

	t_vec planeNorm = m::create<t_vec>({ 0, 0, 1, });
	t_scalar planeD = 5.;

	auto [pos, ty, lam] = m::intersect_line_plane<t_vec>(lineOrg, lineDir, planeNorm, planeD, eps);
	std::cout << "pos: " << pos << ", lam: " << lam << std::endl;
}


template<class t_scalar, class t_vec>
void intersect_plane_poly()
{
	t_scalar eps = 1e-6;

	t_vec planeNorm = m::create<t_vec>({ 1, 0, 0, });
	t_scalar planeD = 0.5;

	std::vector<t_vec> polyVerts
	{{
		m::create<t_vec>({ -1, -1,  5. }),
		m::create<t_vec>({  1, -1,  5. }),
		m::create<t_vec>({  0,  1,  5. }),
	}};

	auto inters = m::intersect_plane_poly<t_vec>(planeNorm, planeD, polyVerts, eps);
	for(const t_vec& pos : inters)
		std::cout << pos << std::endl;
}


int main()
{
	using t_real = double;
	using t_vec = vec<t_real, std::vector>;
	//using t_mat = mat<t_real, std::vector>;

	std::cout << "intersect_line_plane" << std::endl;
	intersect_line_plane<t_real, t_vec>();


	std::cout << "\nintersect_plane_poly" << std::endl;
	intersect_plane_poly<t_real, t_vec>();

	return 0;
}
