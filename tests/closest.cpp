/**
 * closest point test
 * @author Tobias Weber
 * @date 3-aug-2025
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -I../libs -Wall -Wextra -Weffc++ -o closest closest.cpp -lgmp
 */

#include <boost/geometry/index/rtree.hpp>
namespace geo = boost::geometry;
namespace geoidx = geo::index;

#include "matrix_algos.h"
#include "matrix_conts.h"
#include "geo_algos.h"

using namespace m_ops;


int main()
{
	using t_real = double;
	using t_vec = m::vec<t_real, std::vector>;
	using t_mat = m::mat<t_real, std::vector>;

	using t_vertex = geo::model::point<t_real, 3, geo::cs::cartesian>;
	using t_rtree_leaf = std::tuple<t_vertex, std::size_t>;
	using t_rtree = geoidx::rtree<t_rtree_leaf, geoidx::dynamic_linear>;

	const t_real pi = m::pi<t_real>;
	t_mat B = m::B_matrix<t_mat>(4., 5., 6., 90./180.*pi, 90./180.*pi, 60./180.*pi);

	// (hkl) points
	std::vector<t_vec> points;
	points.reserve(8 * 8 * 8);
	for(t_real h = -4.; h < 4.; h += 1.)
	for(t_real k = -4.; k < 4.; k += 1.)
	for(t_real l = -4.; l < 4.; l += 1.)
		points.emplace_back(m::create<t_vec>({ h, k, l }));

	// rtree without and with B
	t_rtree rt = g::make_rtree<
		t_real, t_vec, 3, t_vertex, t_rtree_leaf, t_rtree>(
			points);
	t_rtree rtB = g::make_rtree<
		t_real, t_vec, t_mat, 3, t_vertex, t_rtree_leaf, t_rtree>(
			points, B);

	// query without and with B
	t_vec query_pt = m::create<t_vec>({ 1.5, 1.5, 0. });
	std::vector<std::size_t> pt_indices = g::closest_point_rtree<
		t_real, t_vec, 3, 1, t_vertex, t_rtree_leaf, t_rtree, std::vector>(
			rt, query_pt);
	std::vector<std::size_t> pt_indicesB = g::closest_point_rtree<
		t_real, t_vec, t_mat, 3, 1, t_vertex, t_rtree_leaf, t_rtree, std::vector>(
			rtB, query_pt, B);

	// result without and with B
	std::cout << "direct: ";
	for(std::size_t pt_idx : pt_indices)
		std::cout << points[pt_idx] << std::endl;
	std::cout << "with B: ";
	for(std::size_t pt_idx : pt_indicesB)
		std::cout << points[pt_idx] << std::endl;

	return 0;
}
