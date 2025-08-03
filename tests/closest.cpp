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
	//using t_mat = m::mat<t_real, std::vector>;

	using t_vertex = geo::model::point<t_real, 3, geo::cs::cartesian>;
	using t_rtree_leaf = std::tuple<t_vertex, std::size_t>;
	using t_rtree = geoidx::rtree<t_rtree_leaf, geoidx::dynamic_linear>;

	// (hkl) points
	std::vector<t_vec> points;
	points.reserve(8 * 8 * 8);
	for(t_real h = -4.; h < 4.; h += 1.)
	for(t_real k = -4.; k < 4.; k += 1.)
	for(t_real l = -4.; l < 4.; l += 1.)
		points.emplace_back(m::create<t_vec>({ h, k, l }));

	// rtree
	t_rtree rt = g::make_rtree<
		t_real, t_vec, 3, t_vertex, t_rtree_leaf, t_rtree>(
			points);

	// query
	t_vec query_pt = m::create<t_vec>({ 1., 2., 3 });
	std::vector<std::size_t> pt_indices = g::closest_point_rtree<
		t_real, t_vec, 3, 1, t_vertex, t_rtree_leaf, t_rtree, std::vector>(
			rt, query_pt);

	for(std::size_t pt_idx : pt_indices)
		std::cout << points[pt_idx] << std::endl;

	return 0;
}
