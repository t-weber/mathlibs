/**
 * triangulation tests
 * @author Tobias Weber
 * @date 4-jan-2025
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o triag triag.cpp -lqhullcpp -lqhullstatic_r -lgmp
 */

#include <vector>
#include <iostream>

#include "matrix_algos.h"
#include "geo_algos.h"
#include "matrix_conts.h"

using namespace g;
using namespace m;
using namespace m_ops;

using t_size = std::size_t;
using t_real = double;
using t_vec = vec<t_real, std::vector>;
//using t_mat = mat<t_real, std::vector>;

t_real g_eps = 1e-6;


template<template<class...> class t_cont = std::vector>
void triangulate_2d(const t_cont<t_vec>& vecs)
{
	t_size dim = vecs.begin()->size();

	auto [voronoi1, triags1, neighbours1] = calc_delaunay<t_vec>(dim, vecs, false);
	auto [voronoi2, triags2, neighbours2] = calc_delaunay_iterative<t_vec>(vecs);
	auto [voronoi3, triags3, neighbours3] = calc_delaunay_parabolic<t_vec>(vecs);

	t_size method_idx = 1;
	for(const auto& triags : { triags1, triags2, triags3 })
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "Calculation method " << method_idx << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		// delaunay triangles
		for(const auto& triag : triags)
		{
			std::cout << "Triangle:" << std::endl;
			for(const t_vec& vert : triag)
				std::cout << "\tVertex: " << vert << std::endl;
		}

		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;
		++method_idx;
	}
}


template<template<class...> class t_cont = std::vector>
void triangulate_nd(const t_cont<t_vec>& vecs)
{
	t_size dim = vecs.begin()->size();

	auto [voronoi, polys, neighbours] = calc_delaunay<t_vec>(dim, vecs, false);

		std::cout << "--------------------------------------------------------------------------------" << std::endl;
	// delaunay polygons
	for(const auto& poly : polys)
	{
		std::cout << "Polygon:" << std::endl;
		for(const t_vec& vert : poly)
			std::cout << "\tVertex: " << vert << std::endl;
	}
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
}


int main()
{
	std::vector<t_vec> vecs2d
	{{
		create<t_vec>({ 1.0, 1.1 }),
		create<t_vec>({ 2.0, 1.0 }),
		create<t_vec>({ 2.0, 2.0 }),
		create<t_vec>({ 1.0, 2.1 }),
	}};

	triangulate_2d(vecs2d);


	std::vector<t_vec> vecs3d
	{{
		create<t_vec>({ 1.0, 1.1, 0.5 }),
		create<t_vec>({ 2.0, 1.0, 0.5 }),
		create<t_vec>({ 2.0, 2.0, 0.5 }),
		create<t_vec>({ 1.0, 2.1, 0.5 }),
		create<t_vec>({ 1.5, 5.1, 1.0 }),
	}};

	triangulate_nd(vecs3d);

	return 0;
}
