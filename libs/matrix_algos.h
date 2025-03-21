/**
 * container-agnostic matrix algorithms
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 2017-2023
 * @license see 'LICENSE' file
 *
 * @see references for algorithms:
 * 	- (Arens15): T. Arens et al., ISBN: 978-3-642-44919-2, DOI: 10.1007/978-3-642-44919-2 (2015).
 * 	- (Arfken13): G. B. Arfken et al., ISBN: 978-0-12-384654-9, DOI: 10.1016/C2009-0-30629-7 (2013).
 * 	- (DesktopBronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its HTML version "Desktop Bronstein"].
 * 	- (Bronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its paperback version].
 * 	- (Merziger06): G. Merziger and T. Wirth, ISBN: 3923923333 (2006).
 * 	- (Scarpino11): M. Scarpino, ISBN: 978-1-6172-9017-6 (2011).
 * 	- (Shirane02): G. Shirane et al., ISBN: 978-0-5214-1126-4 (2002).
 * 	- (Kuipers02): J. B. Kuipers, ISBN: 0-691-05872-5 (2002).
 * 	- (FUH 2021): A. Schulz and J. Rollin, "Effiziente Algorithmen", Kurs 1684 (2021), Fernuni Hagen (https://vu.fernuni-hagen.de/lvuweb/lvu/app/Kurs/01684).
 * 	- (Sellers 2014): G. Sellers et al., ISBN: 978-0-321-90294-8 (2014).
 */

#ifndef __MATH_ALGOS_H__
#define __MATH_ALGOS_H__

#include "matrix_concepts.h"

#include <cmath>
#include <complex>
#include <tuple>
#include <vector>
#include <limits>
#include <algorithm>
#include <random>
#include <numeric>
#include <numbers>
#include <optional>
#include <cassert>
#include <iostream>

#ifndef MATH_USE_FLAT_DET
	#define MATH_USE_FLAT_DET 0
#endif

#if MATH_USE_FLAT_DET != 0
	#pragma message("Using flat determinant calculation method.")
#endif


// math
namespace m {

// ----------------------------------------------------------------------------
// forward declarations
// ----------------------------------------------------------------------------
template<class t_vec> requires is_basic_vec<t_vec>
t_vec cross(const t_vec& vec1, const t_vec& vec2);

template<class t_mat, class t_real = typename t_mat::value_type,
	typename t_size = decltype(t_mat{}.size1())>
t_mat givens(t_size N, t_size i, t_size j, t_real angle)
requires is_mat<t_mat>;
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// scalar algos and constants
// ----------------------------------------------------------------------------
template<typename T> constexpr T pi = std::numbers::pi_v<T>;
template<typename T> T golden = std::numbers::phi_v<T>; //T(0.5) + std::sqrt(T(5))/T(2);


/**
 * are two scalars equal within an epsilon range?
 */
template<class T>
bool equals(T t1, T t2, T eps = std::numeric_limits<T>::epsilon())
requires is_scalar<T>
{
	return std::abs(t1 - t2) <= eps;
}


template<typename t_num = unsigned int>
t_num next_multiple(t_num num, t_num granularity)
requires is_scalar<t_num>
{
	t_num div = num / granularity;
	bool rest_is_0 = 1;

	if constexpr(std::is_floating_point_v<t_num>)
	{
		div = std::floor(div);
		t_num rest = std::fmod(num, granularity);
		rest_is_0 = equals(rest, t_num{0});
	}
	else
	{
		t_num rest = num % granularity;
		rest_is_0 = (rest==0);
	}

	return rest_is_0 ? num : (div+1) * granularity;
}


/**
 * mod operation, keeping result positive
 */
template<class t_real>
t_real mod_pos(t_real val, t_real tomod=t_real{2}*pi<t_real>)
requires is_scalar<t_real>
{
	val = std::fmod(val, tomod);
	if(val < t_real(0))
		val += tomod;

	return val;
}


/**
 * are two angles equal within an epsilon range?
 */
template<class T>
bool angle_equals(T t1, T t2, T eps = std::numeric_limits<T>::epsilon(), T tomod=T{2}*pi<T>)
requires is_scalar<T>
{
	t1 = mod_pos<T>(t1, tomod);
	t2 = mod_pos<T>(t2, tomod);

	return std::abs(t1 - t2) <= eps;
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// n-dim algos
// ----------------------------------------------------------------------------
/**
 * are two complex numbers equal within an epsilon range?
 */
template<class T>
bool equals(const T& t1, const T& t2,
	typename T::value_type eps = std::numeric_limits<typename T::value_type>::epsilon())
requires is_complex<T>
{
	return (std::abs(t1.real() - t2.real()) <= eps) &&
		(std::abs(t1.imag() - t2.imag()) <= eps);
}


/**
 * are two vectors equal within an epsilon range?
 */
template<class t_vec>
bool equals(const t_vec& vec1, const t_vec& vec2,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
requires is_basic_vec<t_vec>
{
	using T = typename t_vec::value_type;
	using t_size = decltype(vec1.size());

	// size has to be equal
	if(vec1.size() != vec2.size())
		return false;

	// check each element
	for(t_size i = 0; i < vec1.size(); ++i)
	{
		if constexpr(is_complex<decltype(eps)>)
		{
			if(!equals<T>(vec1[i], vec2[i], eps.real()))
				return false;
		}
		else
		{
			if(!equals<T>(vec1[i], vec2[i], eps))
				return false;
		}
	}

	return true;
}


/**
 * are two matrices equal within an epsilon range?
 */
template<class t_mat>
bool equals(const t_mat& mat1, const t_mat& mat2,
	typename t_mat::value_type eps = std::numeric_limits<typename t_mat::value_type>::epsilon())
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;
	using t_size = decltype(mat1.size1());

	if(mat1.size1() != mat2.size1() || mat1.size2() != mat2.size2())
		return false;

	for(t_size i = 0; i < mat1.size1(); ++i)
	{
		for(t_size j = 0; j < mat1.size2(); ++j)
		{
			if constexpr(is_complex<decltype(eps)>)
			{
				if(!equals<T>(mat1(i,j), mat2(i,j), eps.real()))
					return false;
			}
			else
			{
				if(!equals<T>(mat1(i,j), mat2(i,j), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * are two quaternions equal within an epsilon range?
 * @see (Kuipers02), p. 105
 */
template<class t_quat>
bool equals(const t_quat& quat1, const t_quat& quat2,
	typename t_quat::value_type eps = std::numeric_limits<typename t_quat::value_type>::epsilon())
requires is_basic_quat<t_quat>
{
	using T = typename t_quat::value_type;

	// check each element
	if(!equals<T>(quat1.real(), quat2.real(), eps))
		return false;
	if(!equals<T>(quat1.imag1(), quat2.imag1(), eps))
		return false;
	if(!equals<T>(quat1.imag2(), quat2.imag2(), eps))
		return false;
	if(!equals<T>(quat1.imag3(), quat2.imag3(), eps))
		return false;

	return true;
}


/**
 * create a vector with given size if it is dynamic
 */
template<class t_vec>
t_vec create(decltype(t_vec{}.size()) size = 3)
requires is_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec{size};

	return vec;
}


/**
 * create a matrix with given sizes if it is dynamic
 */
template<class t_mat>
t_mat create(decltype(t_mat{}.size1()) size1, decltype(t_mat{}.size2()) size2)
requires is_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat{size1, size2};

	return mat;
}


/**
 * linearise a matrix to a vector container
 */
template<class t_vec, class t_mat>
t_vec convert(const t_mat& mat)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	//using T_src = typename t_mat::value_type;
	using T_dst = typename t_vec::value_type;
	using t_idx = decltype(mat.size1());

	t_vec vec;

	for(t_idx row_idx = 0; row_idx < mat.size1(); ++row_idx)
		for(t_idx col_idx = 0; col_idx < mat.size2(); ++col_idx)
			vec.push_back(static_cast<T_dst>(mat(row_idx, col_idx)));

	return vec;
}


/**
 * converts a container of objects
 */
template<class t_obj_dst, class t_obj_src, template<class...> class t_cont>
t_cont<t_obj_dst> convert(const t_cont<t_obj_src>& src_objs)
requires (is_vec<t_obj_dst> || is_mat<t_obj_dst>) && (is_vec<t_obj_src> || is_mat<t_obj_src>)
{
	t_cont<t_obj_dst> dst_objs;
	dst_objs.reserve(src_objs.size());

	for(const t_obj_src& src_obj : src_objs)
		dst_objs.emplace_back(convert<t_obj_dst, t_obj_src>(src_obj));

	return dst_objs;
}


/**
 * create a vector from a matrix' diagonal elements
 */
template<class t_vec, class t_mat>
t_vec diagonal(const t_mat& mat)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	//using T_src = typename t_mat::value_type;
	using T_dst = typename t_vec::value_type;
	using t_idx = decltype(mat.size1());

	const t_idx N = std::min(mat.size1(), mat.size2());

	t_vec vec = create<t_vec>(N);
	for(t_idx i = 0; i < N; ++i)
		vec[i] = static_cast<T_dst>(mat(i, i));

	return vec;
}


/**
 * set submatrix to unit
 */
template<class t_mat>
void unit(t_mat& mat,
	decltype(mat.size1()) rows_begin,
	decltype(mat.size2()) cols_begin,
	decltype(mat.size1()) rows_end,
	decltype(mat.size2()) cols_end)
requires is_mat<t_mat>
{
	using t_size = decltype(mat.size1());

	for(t_size i = rows_begin; i < rows_end; ++i)
		for(t_size j = cols_begin; j < cols_end; ++j)
			mat(i,j) = (i==j ? 1 : 0);
}


/**
 * unit matrix
 */
template<class t_mat, typename t_size = decltype(t_mat{}.size1())>
t_mat unit(t_size N1, t_size N2)
requires is_mat<t_mat>
{
	t_mat mat = create<t_mat>(N1, N2);

	unit<t_mat>(mat, 0,0, mat.size1(),mat.size2());
	return mat;
}


/**
 * unit matrix
 */
template<class t_mat, typename t_size = decltype(t_mat{}.size1())>
t_mat unit(t_size N = 0)
requires is_mat<t_mat>
{
	return unit<t_mat>(N,N);
}


/**
 * unit quaternion
 */
template<class t_quat>
const t_quat& unit()
requires is_basic_quat<t_quat>
{
	static const t_quat quat(1, 0, 0, 0);
	return quat;
}


/**
 * converts matrix containers of different value types
 */
template<class t_mat_dst, class t_mat_src>
t_mat_dst convert(const t_mat_src& mat)
requires is_mat<t_mat_dst> && is_mat<t_mat_src>
{
	//using T_src = typename t_mat_src::value_type;
	using T_dst = typename t_mat_dst::value_type;
	using t_idx = decltype(mat.size1());

	t_mat_dst matdst = create<t_mat_dst>(mat.size1(), mat.size2());
	const t_idx act_size1 = std::min(mat.size1(), t_idx(matdst.size1()));
	const t_idx act_size2 = std::min(mat.size2(), t_idx(matdst.size2()));

	for(t_idx row_idx = 0; row_idx < act_size1; ++row_idx)
		for(t_idx col_idx = 0; col_idx < act_size2; ++col_idx)
			matdst(row_idx, col_idx) = static_cast<T_dst>(mat(row_idx, col_idx));

	return matdst;
}


/**
 * converts matrix containers of different value types and possibly sizes
 */
template<class t_mat_dst, class t_mat_src>
void convert(t_mat_dst& mat_dst, const t_mat_src& mat_src)
requires is_mat<t_mat_dst> && is_mat<t_mat_src>
{
	using T_dst = typename t_mat_dst::value_type;
	using t_idx = decltype(mat_src.size1());

	mat_dst = unit<t_mat_dst>(mat_dst.size1(), mat_dst.size2());

	for(t_idx row_idx = 0; row_idx < std::min(mat_src.size1(), mat_dst.size1()); ++row_idx)
		for(t_idx col_idx = 0; col_idx < std::min(mat_src.size2(), mat_dst.size2()); ++col_idx)
			mat_dst(row_idx, col_idx) = static_cast<T_dst>(mat_src(row_idx, col_idx));
}


/**
 * converts vector containers of different value types
 */
template<class t_vec_dst, class t_vec_src>
t_vec_dst convert(const t_vec_src& vec)
requires is_vec<t_vec_dst> && is_vec<t_vec_src>
{
	//using T_src = typename t_vec_src::value_type;
	using T_dst = typename t_vec_dst::value_type;
	using t_idx = decltype(vec.size());

	t_vec_dst vecdst = create<t_vec_dst>(vec.size());

	for(t_idx i = 0; i < vec.size(); ++i)
		vecdst[i] = static_cast<T_dst>(vec[i]);

	return vecdst;
}


/**
 * vector with all values the same
 */
template<class t_vec>
t_vec samevalue(decltype(t_vec{}.size()) N=0, typename t_vec::value_type val=0)
requires is_basic_vec<t_vec>
{
	using size_t = decltype(t_vec{}.size());

	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(N);

	for(size_t i = 0; i < vec.size(); ++i)
		vec[i] = val;

	return vec;
}


/**
 * matrix with all values the same
 */
template<class t_mat, typename t_size = decltype(t_mat{}.size1())>
t_mat samevalue(t_size N1, t_size N2, typename t_mat::value_type val=0)
requires is_mat<t_mat>
{
	t_mat mat = create<t_mat>(N1, N2);

	for(t_size i = 0; i < mat.size1(); ++i)
		for(t_size j = 0; j < mat.size2(); ++j)
			mat(i,j) = val;

	return mat;
}


/**
 * zero matrix
 */
template<class t_mat, typename t_size = decltype(t_mat{}.size1())>
t_mat zero(t_size N1, t_size N2)
requires is_mat<t_mat>
{
	return samevalue<t_mat>(N1, N2, 0);
}


/**
 * zero matrix
 */
template<class t_mat, typename t_size = decltype(t_mat{}.size1())>
t_mat zero(t_size N = 0)
requires is_mat<t_mat>
{
	return zero<t_mat>(N, N);
}


/**
 * zero vector
 */
template<class t_vec>
t_vec zero(decltype(t_vec{}.size()) N=0)
requires is_basic_vec<t_vec>
{
	return samevalue<t_vec>(N, 0);
}



/**
 * zero quaternion
 */
template<class t_quat>
const t_quat& zero()
requires is_basic_quat<t_quat>
{
	static const t_quat quat(0, 0, 0, 0);
	return quat;
}


/**
 * tests for zero scalar
 */
template<class t_scalar>
bool equals_0(t_scalar s, t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
requires is_scalar<t_scalar>
{
	return equals<t_scalar>(s, t_scalar{0}, eps);
}


/**
 * tests for zero complex number
 */
template<class t_cplx, class t_scalar = typename t_cplx::value_type>
bool equals_0(const t_cplx& c, t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
requires is_scalar<t_scalar> && is_complex<t_cplx>
{
	return equals<t_scalar>(c.real(), t_scalar{0}, eps) &&
		equals<t_scalar>(c.imag(), t_scalar{0}, eps);
}


/**
 * tests for zero vector
 */
template<class t_vec>
bool equals_0(const t_vec& vec,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
requires is_basic_vec<t_vec>
{
	return equals<t_vec>(vec, zero<t_vec>(vec.size()), eps);
}


/**
 * tests for zero matrix
 */
template<class t_mat>
bool equals_0(const t_mat& mat,
	typename t_mat::value_type eps = std::numeric_limits<typename t_mat::value_type>::epsilon())
requires is_mat<t_mat>
{
	return equals<t_mat>(mat, zero<t_mat>(mat.size1(), mat.size2()), eps);
}


/**
 * tests for zero quaternion
 */
template<class t_quat>
bool equals_0(const t_quat& quat,
	typename t_quat::value_type eps = std::numeric_limits<typename t_quat::value_type>::epsilon())
requires is_basic_quat<t_quat>
{
	return equals<t_quat>(quat, zero<t_quat>(), eps);
}


/**
 * random scalar
 */
template<class t_scalar>
t_scalar rand(t_scalar min = 0., t_scalar max = 1.)
requires is_scalar<t_scalar>
{
	static std::mt19937 rng{std::random_device{}()};

	if constexpr(std::is_integral_v<t_scalar>)
		return std::uniform_int_distribution<t_scalar>(min, max)(rng);
	else
		return std::uniform_real_distribution<t_scalar>(min, max)(rng);
}


/**
 * random matrix
 */
template<class t_mat, typename t_size = decltype(t_mat{}.size1())>
t_mat rand(t_size N1, t_size N2)
requires is_mat<t_mat>
{
	using t_real = typename t_mat::value_type;

	t_mat mat = create<t_mat>(N1, N2);

	for(t_size i = 0; i < mat.size1(); ++i)
		for(t_size j = 0; j < mat.size2(); ++j)
			mat(i,j) = rand<t_real>();

	return mat;
}


/**
 * random vector
 */
template<class t_vec>
t_vec rand(decltype(t_vec{}.size()) N)
requires is_basic_vec<t_vec>
{
	using t_real = typename t_vec::value_type;
	using t_size = decltype(t_vec{}.size());

	t_vec vec = create<t_vec>(N);

	for(t_size i = 0; i < vec.size(); ++i)
		vec[i] = rand<t_real>();

	return vec;
}


/**
 * random quaternion
 */
template<class t_quat>
const t_quat& rand()
requires is_basic_quat<t_quat>
{
	using t_real = typename t_quat::value_type;

	// TODO: normalise
	return t_quat(
		rand<t_real>(),
		rand<t_real>(),
		rand<t_real>(),
		rand<t_real>());
}


/**
 * tests for symmetric or hermitian matrix
 */
template<class t_mat>
bool is_symm_or_herm(const t_mat& mat,
	typename t_mat::value_type eps = std::numeric_limits<typename t_mat::value_type>::epsilon())
requires is_mat<t_mat>
{
	using t_size = decltype(mat.size1());
	using t_elem = typename t_mat::value_type;

	if(mat.size1() != mat.size2())
		return false;

	for(t_size i = 0; i < mat.size1(); ++i)
	{
		for(t_size j = i + 1; j < mat.size2(); ++j)
		{
			if constexpr(is_complex<t_elem>)
			{
				// not hermitian?
				if(!equals<t_elem>(mat(i,j), std::conj(mat(j,i)), eps))
					return false;
			}
			else
			{
				// not symmetric?
				if(!equals<t_elem>(mat(i,j), mat(j,i), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * transpose matrix
 * WARNING: not possible for static non-square matrix!
 */
template<class t_mat>
t_mat trans(const t_mat& mat)
requires is_mat<t_mat>
{
	using t_size = decltype(mat.size1());

	t_mat mat2 = create<t_mat>(mat.size2(), mat.size1());

	for(t_size i = 0; i < mat.size1(); ++i)
		for(t_size j = 0; j < mat.size2(); ++j)
			mat2(j,i) = mat(i,j);

	return mat2;
}


/**
 * create vector from initializer_list
 */
template<class t_vec>
t_vec create(const std::initializer_list<typename t_vec::value_type>& lst)
requires is_basic_vec<t_vec>
{
	t_vec vec = create<t_vec>(lst.size());
	using t_size = decltype(vec.size());

	auto iterLst = lst.begin();
	for(t_size i = 0; i < vec.size(); ++i)
	{
		if(iterLst != lst.end())
		{
			vec[i] = *iterLst;
			std::advance(iterLst, 1);
		}
		else	// vector larger than given list?
		{
			vec[i] = 0;
		}
	}

	return vec;
}


/**
 * create matrix from nested initializer_lists in columns/rows order
 */
template<class t_mat,
	template<class...> class t_cont_outer = std::initializer_list,
	template<class...> class t_cont = std::initializer_list>
t_mat create(const t_cont_outer<t_cont<typename t_mat::value_type>>& lst)
requires is_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());

	const t_size col_idxs = lst.size();
	const t_size row_idxs = lst.begin()->size();

	t_mat mat = unit<t_mat>(row_idxs, col_idxs);

	auto iterCol = lst.begin();
	for(t_size col_idx = 0; col_idx < col_idxs; ++col_idx)
	{
		auto iterRow = iterCol->begin();
		for(t_size row_idx = 0; row_idx < row_idxs; ++row_idx)
		{
			mat(row_idx, col_idx) = *iterRow;
			std::advance(iterRow, 1);
		}

		std::advance(iterCol, 1);
	}

	return mat;
}


/**
 * create matrix from column (or row) vectors
 */
template<class t_mat, class t_vec,
	template<class...> class t_cont_outer = std::initializer_list>
t_mat create(const t_cont_outer<t_vec>& lst, bool bRow = false)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	using t_size = decltype(t_mat{}.size1());

	const t_size col_idxs = lst.size();
	const t_size row_idxs = lst.begin()->size();

	t_mat mat = unit<t_mat>(row_idxs, col_idxs);

	auto iterCol = lst.begin();
	for(t_size col_idx = 0; col_idx < col_idxs; ++col_idx)
	{
		for(t_size row_idx = 0; row_idx < row_idxs; ++row_idx)
			mat(row_idx, col_idx) = (*iterCol)[row_idx];
		std::advance(iterCol, 1);
	}

	if(bRow) mat = trans<t_mat>(mat);
	return mat;
}


/**
 * create matrix from initializer_list in column/row order
 */
template<class t_mat>
t_mat create(const std::initializer_list<typename t_mat::value_type>& lst)
requires is_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());
	const t_size N = std::sqrt(lst.size());

	t_mat mat = unit<t_mat>(N, N);

	auto iter = lst.begin();
	for(t_size row_idx = 0; row_idx < N; ++row_idx)
	{
		for(t_size col_idx = 0; col_idx < N; ++col_idx)
		{
			mat(row_idx, col_idx) = *iter;
			std::advance(iter, 1);
		}
	}

	return mat;
}


/**
 * create matrix from initializer_list in column/row order
 */
template<class t_mat, class t_init =  std::initializer_list<typename t_mat::value_type>,
	class t_size = decltype(t_mat{}.size1())>
t_mat create(t_size rows, t_size cols, const t_init& lst)
requires is_mat<t_mat>
{
	t_mat mat = zero<t_mat>(rows, cols);

	auto iter = lst.begin();
	for(t_size row_idx = 0; row_idx < rows; ++row_idx)
	{
		for(t_size col_idx = 0; col_idx < cols; ++col_idx)
		{
			mat(row_idx, col_idx) = *iter;
			std::advance(iter, 1);
		}
	}

	return mat;
}


/**
 * get a column vector from a matrix
 */
template<class t_mat, class t_vec, typename t_size = decltype(t_mat{}.size1())>
t_vec col(const t_mat& mat, t_size col)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec = create<t_vec>(mat.size1());

	for(t_size i = 0; i < static_cast<t_size>(mat.size1()); ++i)
		vec[i] = mat(i, col);

	return vec;
}


/**
 * get a row vector from a matrix
 */
template<class t_mat, class t_vec, typename t_size = decltype(t_mat{}.size1())>
t_vec row(const t_mat& mat, t_size row)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec = create<t_vec>(mat.size2());

	for(t_size i = 0; i < mat.size2(); ++i)
		vec[i] = mat(row, i);

	return vec;
}


/**
 * set a column vector in a matrix
 */
template<class t_mat, class t_vec, typename t_size = decltype(t_mat{}.size1())>
void set_col(t_mat& mat, const t_vec& vec, t_size col)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	t_size num_rows = std::min<t_size>(mat.size1(), static_cast<t_size>(vec.size()));
	for(t_size i = 0; i < num_rows; ++i)
		mat(i, col) = vec[i];
}


/**
 * set a row vector in a matrix
 */
template<class t_mat, class t_vec, typename t_size = decltype(t_mat{}.size1())>
void set_row(t_mat& mat, const t_vec& vec, t_size row)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	for(t_size i = 0; i < std::min(mat.size2(), vec.size()); ++i)
		mat(row, i) = vec[i];
}


/**
 * inner product <vec1|vec2>
 */
template<class t_vec>
typename t_vec::value_type inner(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>
{
	using t_size = decltype(vec1.size());
	typename t_vec::value_type val(0);

	const t_size I = vec1.size();
	[[__maybe_unused__]] const t_size J = vec2.size();
	assert(I == J);

	for(t_size i = 0; i < I; ++i)
	{
		if constexpr(is_complex<typename t_vec::value_type>)
			val += std::conj(vec1[i]) * vec2[i];
		else
			val += vec1[i] * vec2[i];
	}

	return val;
}


/**
 * inner product between two vectors of different type
 */
template<class t_vec1, class t_vec2>
typename t_vec1::value_type inner(const t_vec1& vec1, const t_vec2& vec2)
requires is_basic_vec<t_vec1> && is_basic_vec<t_vec2>
{
	using t_size = decltype(vec1.size());

	const t_size I = vec1.size();
	const t_size J = vec2.size();
	//assert(I == J);

	if(I==0 || J==0)
		return typename t_vec1::value_type{};

	// first element
	auto val = vec1[0]*vec2[0];

	// remaining elements
	for(t_size i = 1; i < std::min(I, J); ++i)
	{
		if constexpr(is_complex<typename t_vec1::value_type>)
		{
			auto prod = std::conj(vec1[i]) * vec2[i];
			val = val + prod;
		}
		else
		{
			auto prod = vec1[i]*vec2[i];
			val = val + prod;
		}
	}

	return val;
}


/**
 * sum components of a vector
 */
template<class t_vec>
typename t_vec::value_type sum(const t_vec& vec)
requires is_basic_vec<t_vec>
{
	using t_size = decltype(vec.size());
	typename t_vec::value_type val(0);

	for(t_size i = 0; i < vec.size(); ++i)
		val += vec[i];

	return val;
}


/**
 * 2-norm
 */
template<class t_vec>
typename t_vec::value_type norm(const t_vec& vec)
requires is_basic_vec<t_vec>
{
	return std::sqrt(inner<t_vec>(vec, vec));
}


/**
 * n-norm
 */
template<class t_vec, class t_real = typename t_vec::value_type>
typename t_vec::value_type norm(const t_vec& vec, t_real n)
requires is_basic_vec<t_vec>
{
	using t_size = decltype(vec.size());

	t_real d = t_real{0};
	for(t_size i = 0; i < vec.size(); ++i)
		d += std::pow(std::abs(vec[i]), n);
	n = std::pow(d, t_real(1)/n);
	return n;
}


/**
 * outer product |v1><v2|
 */
template<class t_mat, class t_vec>
t_mat outer(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(vec1.size());

	const t_size N1 = vec1.size();
	const t_size N2 = vec2.size();
	t_mat mat = create<t_mat>(N1, N2);

	for(t_size n1 = 0; n1 < N1; ++n1)
	{
		for(t_size n2 = 0; n2 < N2; ++n2)
		{
			if constexpr(is_complex<typename t_vec::value_type>)
				mat(n1, n2) = std::conj(vec1[n1]) * vec2[n2];
			else
				mat(n1, n2) = vec1[n1]*vec2[n2];
		}
	}

	return mat;
}


/**
 * outer product |v1><v2|, "flattened" to a (state) vector
 */
template<class t_vec, class t_mat>
t_vec outer_flat(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(vec1.size());
	const t_size ROWS = vec1.size();
	const t_size COLS = vec2.size();

	t_mat outer = m::outer<t_mat, t_vec>(vec1, vec2);
	t_vec outer_flat = create<t_vec>(ROWS*COLS);

	for(t_size i = 0; i < ROWS; ++i)
		for(t_size j = 0; j < COLS; ++j)
			outer_flat[i*COLS + j] = outer(i, j);

	return outer_flat;
}


/**
 * outer products of a collection of vectors
 */
template<class t_vec, class t_mat, template<class...> class t_cont = std::initializer_list>
t_vec outer_flat(const t_cont<t_vec>& vecs)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	if(vecs.size() == 0)
		return t_vec{};
	else if(vecs.size() == 1)
		return *vecs.begin();

	auto iter = vecs.begin();
	t_vec outer = outer_flat<t_vec, t_mat>(*iter, *std::next(iter, 1));

	for(iter = std::next(iter, 2); iter != vecs.end(); std::advance(iter, 1))
		outer = outer_flat<t_vec, t_mat>(outer, *iter);

	return outer;
}


/**
 * outer/tensor product
 */
template<class t_mat>
t_mat outer(const t_mat& mat1, const t_mat& mat2)
requires is_mat<t_mat>
{
	using t_size = decltype(mat1.size1());

	const t_size m1s1 = mat1.size1();
	const t_size m1s2 = mat1.size2();
	const t_size m2s1 = mat2.size1();
	const t_size m2s2 = mat2.size2();
	t_mat mat = create<t_mat>(m1s1*m2s1, m1s2*m2s2);

	for(t_size i1 = 0; i1 < m1s1; ++i1)
		for(t_size j1 = 0; j1 < m1s2; ++j1)
			for(t_size i2 = 0; i2 < m2s1; ++i2)
				for(t_size j2 = 0; j2 < m2s2; ++j2)
					mat(i1*m2s1+i2, j1*m2s2+j2) = mat1(i1, j1) * mat2(i2, j2);

	return mat;
}


/**
 * matrix-matrix product, M_ij = R_ik S_kj
 */
template<class t_mat>
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
t_mat mult(const t_mat& R, const t_mat& S) noexcept
{
	if constexpr(m::is_dyn_mat<t_mat>)
		assert((R.size2() == S.size1()));
	else
		static_assert(R.size2() == S.size1());

	//using t_scalar = std::common_type_t<t_scalar_1, t_scalar_2>;
	using t_scalar = typename t_mat::value_type;
	using t_size = decltype(t_mat{}.size1());

	const t_size I = R.size1();
	const t_size J = S.size2();
	const t_size K = R.size2();

	t_mat M = create<t_mat>(I, J);

	for(t_size i = 0; i < I; ++i)
	{
		for(t_size j = 0; j < J; ++j)
		{
			M(i, j) = t_scalar{0};

			for(t_size k = 0; k < K; ++k)
				M(i, j) += R(i, k) * S(k, j);
		}
	}

	return M;
}


/**
 * matrix-vector product
 */
template<class t_mat, class t_vec>
t_vec mult(const t_mat& mat, const t_vec& vec)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
	&& m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	using t_size = decltype(t_mat{}.size1());

	if constexpr(m::is_dyn_mat<t_mat>)
		assert((mat.size2() == t_size(vec.size())));
	else
		static_assert(mat.size2() == t_size(vec.size()));

	t_vec vecRet = m::create<t_vec>(mat.size1());

	for(t_size row = 0; row < mat.size1(); ++row)
	{
		vecRet[row] = typename t_vec::value_type{/*0*/};
		for(t_size col = 0; col < mat.size2(); ++col)
		{
			auto elem = mat(row, col) * vec[col];
			vecRet[row] += elem;
		}
	}

	return vecRet;
}


/**
 * matrix-vector product using only a portion of the matrix
 */
template<class t_mat, class t_vec, typename t_size = decltype(t_mat{}.size1())>
t_vec mult(const t_mat& mat, const t_vec& vec,
	t_size outsize /*= std::numeric_limits<t_size>::max()*/,
	t_size row_begin = 0, t_size col_begin = 0)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
	&& m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	t_size insize = std::min(vec.size(), mat.size2()-col_begin);
	outsize = std::min(outsize, mat.size1()-row_begin);
	t_vec vecRet(outsize);

	for(t_size row = row_begin; row < row_begin+outsize; ++row)
	{
		vecRet[row-row_begin] = t_real{/*0*/};
		for(t_size col = col_begin; col < col_begin+insize; ++col)
			vecRet[row-row_begin] += mat(row, col) * vec[col-col_begin];
	}

	return vecRet;
}


/**
 * component-wise min/max, e.g. for bounding box
 */
template<class t_vec, template<class...> class t_cont = std::initializer_list>
std::pair<t_vec, t_vec> minmax_comp(const t_cont<t_vec>& vecs)
requires is_basic_vec<t_vec>
{
	if(vecs.size() == 0)
		return std::make_pair(t_vec{}, t_vec{});
	else if(vecs.size() == 1)
		return std::make_pair(*vecs.begin(), *vecs.begin());

	using t_size = decltype(vecs[0].size());

	const t_size dim = vecs[0].size();
	t_vec min_vec = create<t_vec>(dim);
	t_vec max_vec = create<t_vec>(dim);

	for(t_size comp_idx = 0; comp_idx < dim; ++comp_idx)
	{
		auto [min_iter, max_iter] = std::minmax_element(vecs.begin(), vecs.end(),
			[comp_idx](const t_vec& vec0, const t_vec& vec1) -> bool
		{
			return vec0[comp_idx] < vec1[comp_idx];
		});

		min_vec[comp_idx] = (*min_iter)[comp_idx];
		max_vec[comp_idx] = (*max_iter)[comp_idx];
	}

	return std::make_pair(min_vec, max_vec);
}


/**
 * get min/max distance from a given point
 */
template<class t_vec, template<class...> class t_cont = std::initializer_list>
std::pair<typename t_vec::value_type, typename t_vec::value_type>
minmax_dist(const t_cont<t_vec>& vecs, const t_vec& pt)
requires is_basic_vec<t_vec>
{
	if(vecs.size() == 0)
		return std::make_pair(0, 0);

	using t_real = typename t_vec::value_type;

	std::vector<t_real> dists(vecs.size());
	std::transform(vecs.begin(), vecs.end(), dists.begin(),
		[&pt](const t_vec& vec) -> t_real
	{
		t_vec diff_vec = vec - pt;
		return inner<t_vec>(diff_vec, diff_vec);
	});

	auto [min_dist, max_dist] = std::minmax_element(dists.begin(), dists.end());
	return std::make_pair(std::sqrt(*min_dist), std::sqrt(*max_dist));
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// with metric
// ----------------------------------------------------------------------------

/**
 * covariant metric tensor, g_{i,j} = e_i * e_j
 * @see (Arens15), p. 808
 */
template<class t_mat, class t_vec, template<class...> class t_cont=std::initializer_list>
t_mat metric(const t_cont<t_vec>& basis_co)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	using t_size = decltype(t_mat{}.size1());

	const t_size N = basis_co.size();
	t_mat g_co = create<t_mat>(N, N);

	auto iter_i = basis_co.begin();
	for(t_size i = 0; i < N; ++i)
	{
		auto iter_j = basis_co.begin();
		for(t_size j = 0; j < N; ++j)
		{
			g_co(i,j) = inner<t_vec>(*iter_i, *iter_j);
			std::advance(iter_j, 1);
		}
		std::advance(iter_i, 1);
	}

	return g_co;
}


/**
 * lower index using metric
 * @see (Arens15), p. 808
 */
template<class t_mat, class t_vec>
t_vec lower_index(const t_mat& metric_co, const t_vec& vec_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	using t_size = decltype(vec_contra.size());

	const t_size N = vec_contra.size();
	t_vec vec_co = zero<t_vec>(N);

	for(t_size i = 0; i < N; ++i)
		for(t_size j = 0; j < N; ++j)
			vec_co[i] += metric_co(i,j) * vec_contra[j];

	return vec_co;
}


/**
 * raise index using metric
 * @see (Arens15), p. 808
 */
template<class t_mat, class t_vec>
t_vec raise_index(const t_mat& metric_contra, const t_vec& vec_co)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	using t_size = decltype(vec_co.size());

	const t_size N = vec_co.size();
	t_vec vec_contra = zero<t_vec>(N);

	for(t_size i = 0; i < N; ++i)
		for(t_size j = 0; j < N; ++j)
			vec_contra[i] += metric_contra(i,j) * vec_co[j];

	return vec_contra;
}


/**
 * inner product using metric
 * @see (Arens15), p. 808
 */
template<class t_mat, class t_vec>
typename t_vec::value_type inner(const t_mat& metric_co, const t_vec& vec1_contra, const t_vec& vec2_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec2_co = lower_index<t_mat, t_vec>(metric_co, vec2_contra);
	return inner<t_vec>(vec1_contra, vec2_co);
}


/**
 * 2-norm using metric
 * @see (Arens15), p. 808
 */
template<class t_mat, class t_vec>
typename t_vec::value_type norm(const t_mat& metric_co, const t_vec& vec_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	return std::sqrt(inner<t_mat, t_vec>(metric_co, vec_contra, vec_contra));
}
// ----------------------------------------------------------------------------




/**
 * matrix to project onto vector: P = |v><v|
 * from: |x'> = <v|x> * |v> = |v><v|x> = |v><v| * |x>
 * @see (Arens15), p. 814
 */
template<class t_mat, class t_vec>
t_mat projector(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	if(is_normalised)
	{
		return outer<t_mat, t_vec>(vec, vec);
	}
	else
	{
		const auto len = norm<t_vec>(vec);
		t_vec _vec = vec / len;
		return outer<t_mat, t_vec>(_vec, _vec);
	}
}


/**
 * project vector vec onto another vector vecProj
 * @see (Arens15), p. 814
 */
template<class t_vec>
t_vec project(const t_vec& vec, const t_vec& vecProj, bool is_normalised = true)
requires is_vec<t_vec>
{
	if(is_normalised)
	{
		return inner<t_vec>(vecProj, vec) * vecProj;
	}
	else
	{
		const auto len = norm<t_vec>(vecProj);
		const t_vec _vecProj = vecProj / len;
		return inner<t_vec>(_vecProj, vec) * _vecProj;
	}
}


/**
 * project vector vec onto another vector vecProj
 * don't multiply with direction vector
 * @see (Arens15), p. 814
 */
template<class t_vec>
typename t_vec::value_type
project_scalar(const t_vec& vec, const t_vec& vecProj, bool is_normalised = true)
requires is_vec<t_vec>
{
	if(is_normalised)
	{
		return inner<t_vec>(vecProj, vec);
	}
	else
	{
		const auto len = norm<t_vec>(vecProj);
		const t_vec _vecProj = vecProj / len;
		return inner<t_vec>(_vecProj, vec);
	}
}


/**
 * project vector vec onto the line lineOrigin + lam*lineDir
 * shift line to go through origin, calculate projection and shift back
 * @returns [closest point, distance]
 */
template<class t_vec>
std::tuple<t_vec, typename t_vec::value_type>
project_line(const t_vec& vec,
	const t_vec& lineOrigin, const t_vec& lineDir, bool is_normalised = true)
requires is_vec<t_vec>
{
	const t_vec ptShifted = vec - lineOrigin;
	const t_vec ptProj = project<t_vec>(ptShifted, lineDir, is_normalised);
	const t_vec ptNearest = lineOrigin + ptProj;

	const typename t_vec::value_type dist = norm<t_vec>(vec - ptNearest);
	return std::make_tuple(ptNearest, dist);
}


/**
 * distance between point and line
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_real dist_pt_line(const t_vec& pt,
	const t_vec& linePt1, const t_vec& linePt2,
	bool bLineIsFinite=true)
requires is_vec<t_vec>
{
	using t_size = decltype(pt.size());
	const t_size dim = linePt1.size();

	const t_vec lineDir = linePt2 - linePt1;
	const auto [nearestPt, dist] = project_line<t_vec>(pt, linePt1, lineDir, false);


	// get point component with max. difference
	t_real diff = -1.;
	t_size compidx = 0;
	for(t_size i = 0; i < dim; ++i)
	{
		t_real newdiff = std::abs(linePt2[i] - linePt1[i]);
		if(newdiff > diff)
		{
			diff = newdiff;
			compidx = i;
		}
	}


	t_real t = (nearestPt[compidx]-linePt1[compidx]) / (linePt2[compidx]-linePt1[compidx]);
	if(bLineIsFinite && t>=t_real{0} && t<=t_real{1})
	{
		// projection is on line -> use distance between point and projection
		return dist;
	}
	else
	{
		// projection is not on line -> use distance between point and closest line end point
		if(std::abs(t-t_real{0}) < std::abs(t-t_real{1}))
			return norm<t_vec>(linePt1 - pt);
		else
			return norm<t_vec>(linePt2 - pt);
	}
}


/**
 * matrix to project onto orthogonal complement (plane perpendicular to vector): P = 1-|v><v|
 * from completeness relation: 1 = sum_i |v_i><v_i| = |x><x| + |y><y| + |z><z|
 * @see (Arens15), p. 814
 */
template<class t_mat, class t_vec>
t_mat ortho_projector(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(vec.size());

	const t_size iSize = vec.size();
	return unit<t_mat>(iSize) -
		projector<t_mat, t_vec>(vec, is_normalised);
}


/**
 * matrix to mirror on plane perpendicular to vector: P = 1 - 2*|v><v|
 * subtracts twice its projection onto the plane normal from the vector
 * @see (Arens15), p. 710
 *
 * this operation is used for the grover iterations
 * @see (FUH 2021), p. 26f.
 */
template<class t_mat, class t_vec>
t_mat ortho_mirror_op(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(vec.size());
	using T = typename t_vec::value_type;

	const t_size iSize = vec.size();

	return unit<t_mat>(iSize) -
		T(2)*projector<t_mat, t_vec>(vec, is_normalised);
}


/**
 * matrix to mirror [a, b, c, ...] into, e.g.,  [a, b', 0, 0]
 * @see (Scarpino11), p. 268
 */
template<class t_mat, class t_vec>
std::tuple<t_mat, bool> ortho_mirror_zero_op(const t_vec& vec, decltype(vec.size()) row)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(vec.size());
	using T = typename t_vec::value_type;

	const t_size N = vec.size();

	t_vec vecSub = zero<t_vec>(N);
	for(t_size i = 0; i < row; ++i)
		vecSub[i] = vec[i];

	// norm of rest vector
	T n = T(0);
	for(t_size i = row; i < N; ++i)
		n += vec[i]*vec[i];
	vecSub[row] = std::sqrt(n);

	const t_vec vecOp = vec - vecSub;

	// nothing to do -> return unit matrix
	if(equals_0<t_vec>(vecOp))
		return std::make_tuple(unit<t_mat>(vecOp.size(), vecOp.size()), false);

	return std::make_tuple(ortho_mirror_op<t_mat, t_vec>(vecOp, false), true);
}


/**
 * QR decomposition of a matrix
 * @returns [Q, R, number of mirror operations]
 * @see (Scarpino11), pp. 269-272
 */
template<class t_mat, class t_vec>
std::tuple<t_mat, t_mat, decltype(t_mat{}.size1())> qr(const t_mat& mat)
requires is_mat<t_mat> && is_vec<t_vec>
{
	//using T = typename t_mat::value_type;
	using size_t = decltype(mat.size1());

	const size_t rows = mat.size1();
	const size_t cols = mat.size2();
	const size_t N = std::min(cols, rows);

	t_mat R = mat;
	t_mat Q = unit<t_mat>(N, N);

	size_t num_mirrors = 0;
	for(size_t icol = 0; icol < N-1; ++icol)
	{
		t_vec vecCol = col<t_mat, t_vec>(R, icol);
		const auto [matMirror, reflected] = ortho_mirror_zero_op<t_mat, t_vec>(vecCol, icol);
		Q = Q * matMirror;
		R = matMirror * R;

		if(reflected)
			++num_mirrors;
	}

	return std::make_tuple(Q, R, num_mirrors);
}


/**
 * project vector vec onto plane through the origin and perpendicular to vector vecNorm
 * (e.g. used to calculate magnetic interaction vector M_perp)
 */
template<class t_vec>
t_vec ortho_project(const t_vec& vec, const t_vec& vecNorm, bool is_normalised = true)
requires is_vec<t_vec>
{
	//const std::size_t iSize = vec.size();
	return vec - project<t_vec>(vec, vecNorm, is_normalised);
}


/**
 * project vector vec onto plane perpendicular to vector vecNorm with distance d
 * vecNorm has to be normalised and plane in Hessian form: x*vecNorm = d
 */
template<class t_vec>
t_vec ortho_project_plane(const t_vec& vec,
	const t_vec& vecNorm, typename t_vec::value_type d)
requires is_vec<t_vec>
{
	// project onto plane through origin
	t_vec vecProj0 = ortho_project<t_vec>(vec, vecNorm, 1);
	// add distance of plane to origin
	return vecProj0 + d*vecNorm;
}


/**
 * mirror a vector on a plane perpendicular to vector vecNorm with distance d
 * vecNorm has to be normalised and plane in Hessian form: x*vecNorm = d
 * @see (Arens15), p. 710
 */
template<class t_vec>
t_vec ortho_mirror_plane(const t_vec& vec,
	const t_vec& vecNorm, typename t_vec::value_type d)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vecProj = ortho_project_plane<t_vec>(vec, vecNorm, d);
	return vec - T(2)*(vec - vecProj);
}


/**
 * find orthonormal substitute basis for vector space (Gram-Schmidt algo)
 * remove orthogonal projections to all other basis vectors: |i'> = (1 - sum_{j<i} |j><j|) |i>
 * @see (Arens15), p. 744
 * @see https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
 */
template<class t_vec,
	template<class...> class t_cont_out = std::vector,
	template<class...> class t_cont_in = std::initializer_list>
t_cont_out<t_vec> orthonorm_sys(const t_cont_in<t_vec>& sys)
requires is_vec<t_vec>
{
	t_cont_out<t_vec> newsys;
	newsys.reserve(sys.size());

	//const std::size_t N = sys.size();
	for(const t_vec& vecSys : sys)
	{
		t_vec vecOrthoProj = vecSys;

		// subtract projections to other basis vectors
		for(const t_vec& vecNewSys : newsys)
			vecOrthoProj -= project<t_vec>(vecSys, vecNewSys, true);

		// normalise
		vecOrthoProj /= norm<t_vec>(vecOrthoProj);
		newsys.emplace_back(std::move(vecOrthoProj));
	}

	return newsys;
}


/**
 * find orthonormal substitute basis for vector space (Gram-Schmidt algo)
 * remove orthogonal projections to all other basis vectors: |i'> = (1 - sum_{j<i} |j><j|) |i>
 * @see (Arens15), p. 744
 * @see https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
 */
template<class t_mat, class t_vec>
t_mat orthonorm(const t_mat& mat)
requires is_mat<t_mat> && is_vec<t_vec>
{
	using t_size = decltype(mat.size1());
	//using t_real = typename t_mat::value_type;

	t_mat matOut = mat;

	for(t_size colidx = 0; colidx < mat.size2(); ++colidx)
	{
		t_vec vecSys = col<t_mat, t_vec>(mat, colidx);
		t_vec vecOrthoProj = vecSys;

		// subtract projections to other basis vectors
		for(t_size newcolidx = 0; newcolidx < colidx; ++newcolidx)
			vecOrthoProj -= project<t_vec>(vecSys, col<t_mat, t_vec>(matOut, newcolidx), true);

		// normalise and set column
		vecOrthoProj /= norm<t_vec>(vecOrthoProj);
		set_col(matOut, vecOrthoProj, colidx);
	}

	return matOut;
}


/**
 * submatrix removing a column/row from a matrix stored in a vector container
 */
template<class t_vec, class t_matvec = t_vec, typename t_size = decltype(t_vec{}.size())>
t_vec flat_submat(const t_matvec& mat,
	t_size iNumRows, t_size iNumCols,
	t_size rem_row, t_size rem_col)
requires is_basic_vec<t_vec>
{
	t_vec vec;
	vec.reserve(iNumRows);

	for(t_size row_idx = 0; row_idx < iNumRows; ++row_idx)
	{
		if(row_idx == rem_row)
			continue;

		for(t_size col_idx = 0; col_idx < iNumCols; ++col_idx)
		{
			if(col_idx == rem_col)
				continue;
			vec.push_back(mat[row_idx*iNumCols + col_idx]);
		}
	}

	return vec;
}


/**
 * submatrix removing a column/row from a matrix
 */
template<class t_mat>
t_mat submat(const t_mat& mat, decltype(mat.size1()) rem_row, decltype(mat.size2()) rem_col)
requires is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());
	t_mat matRet = m::create<t_mat>(mat.size1()-1, mat.size2()-1);

	size_t res_row = 0;
	for(size_t row_idx = 0; row_idx < mat.size1(); ++row_idx)
	{
		if(row_idx == rem_row)
			continue;

		size_t iResCol = 0;
		for(size_t col_idx = 0; col_idx < mat.size2(); ++col_idx)
		{
			if(col_idx == rem_col)
				continue;

			matRet(res_row, iResCol) = mat(row_idx, col_idx);
			++iResCol;
		}

		++res_row;
	}

	return matRet;
}


/**
 * determinant from a square matrix stored in a vector container
 * @see (Merziger06), p. 185
 */
template<class t_vec, class t_matvec = t_vec, typename t_size = decltype(t_vec{}.size())>
typename t_vec::value_type flat_det(const t_matvec& mat, t_size N)
requires is_basic_vec<t_vec>
{
	using T = typename t_vec::value_type;

	// special cases
	if(N == 0)
		return 0;
	else if(N == 1)
		return mat[0];
	else if(N == 2)
		return mat[0]*mat[3] - mat[1]*mat[2];


	T fullDet = T(0);
	t_size row_idx = 0;

	// get row with maximum number of zeros
	t_size max_num_zeros = 0;
	for(t_size cur_row = 0; cur_row < N; ++cur_row)
	{
		t_size num_zeros = 0;
		for(t_size cur_col = 0; cur_col < N; ++cur_col)
		{
			if(equals<T>(mat[cur_row*N + cur_col], T(0)))
				++num_zeros;
		}

		if(num_zeros > max_num_zeros)
		{
			row_idx = cur_row;
			max_num_zeros = num_zeros;
		}
	}


	// recursively expand determinant along a row
	for(t_size col_idx = 0; col_idx < N; ++col_idx)
	{
		const T elem = mat[row_idx*N + col_idx];
		if(equals<T>(elem, 0))
			continue;

		const T sgn = ((row_idx+col_idx) % 2) == 0 ? T(1) : T(-1);
		const t_vec subMat = flat_submat<t_vec, t_matvec>(mat, N, N, row_idx, col_idx);
		const T subDet = flat_det<t_vec>(subMat, N-1) * sgn;

		fullDet += elem * subDet;
	}

	return fullDet;
}


/**
 * determinant
 * @see (Merziger06), p. 185
 */
template<class t_mat, class t_vec>
typename t_mat::value_type det(const t_mat& mat)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;
	using size_t = decltype(mat.size1());

	if(mat.size1() != mat.size2())
		return 0;

	size_t N = mat.size1();

	// special cases
	switch(N)
	{
		case 0: return 0;
		case 1: return mat(0, 0);
		case 2: return mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);
	}

	T res = T{1};

#if MATH_USE_FLAT_DET == 0
	const auto [Q, R, num_mirrors] = qr<t_mat, t_vec>(mat);

	for(size_t i = 0; i < N; ++i)
		res *= R(i,i);

	// odd number of mirror operations for qr
	if((num_mirrors % 2) != 0)
		res = -res;

	// test sign of det(Q)
	//std::vector<T> matFlatQ = convert<std::vector<T>, t_mat>(Q);
	//T detQ = flat_det<std::vector<T>>(matFlatQ, Q.size1());
	//if(detQ < 0.) res = -res;

#else
	//std::vector<T> matFlat = convert<std::vector<T>, t_mat>(mat);
	const auto& matFlat = matvec_adapter<t_mat>{mat};
	res = flat_det<t_vec>(matFlat, mat.size1());
#endif

	return res;
}


/**
 * trace
 */
template<class t_mat>
typename t_mat::value_type trace(const t_mat& mat)
requires is_mat<t_mat>
{
	using t_size = decltype(mat.size1());
	using T = typename t_mat::value_type;

	T _tr = T(0);

	t_size N = std::min(mat.size1(), mat.size2());
	for(t_size i = 0; i < N; ++i)
		_tr += mat(i,i);

	return _tr;
}


/**
 * inverted matrix
 * @see https://en.wikipedia.org/wiki/Invertible_matrix#In_relation_to_its_adjugate
 * @see https://en.wikipedia.org/wiki/Adjugate_matrix
 */
template<class t_mat, class t_vec>
std::tuple<t_mat, bool> inv(const t_mat& mat)
requires is_mat<t_mat> && is_vec<t_vec>
{
	using T = typename t_mat::value_type;
	using t_idx = decltype(mat.size1());

	const t_idx N = mat.size1();

	// fail if matrix is not square
	if(N != mat.size2())
		return std::make_tuple(t_mat(), false);

#if MATH_USE_FLAT_DET == 0
	const T fullDet = det<t_mat, t_vec>(mat);

#else
	//using t_matvec = std::vector<T>;
	//const t_matvec matFlat = convert<t_matvec, t_mat>(mat);
	const auto& matFlat = matvec_adapter<t_mat>{mat};
	const T fullDet = flat_det<t_vec>(matFlat, N);
#endif

	// fail if determinant is zero
	if(equals<T>(fullDet, 0))
		return std::make_tuple(t_mat(), false);

	t_mat matInv = create<t_mat>(N, N);

	for(t_idx i = 0; i < N; ++i)
	{
		for(t_idx j = 0; j < N; ++j)
		{
#if MATH_USE_FLAT_DET == 0
			// careful with size1() and size2() of static matrices (not fulfilling is_dyn_mat!
			const t_mat subMat = submat<t_mat>(mat, i, j);
			const T subDet = det<t_mat, t_vec>(subMat);

#else
			// alternatively, better for static matrices:
			const t_vec subMat = flat_submat<t_vec>(matFlat, N, N, i, j);
			const T subDet = flat_det<t_vec>(subMat, N-1);
#endif

			const T sgn = ((i+j) % 2) == 0 ? T(1) : T(-1);
			matInv(j,i) = sgn * subDet;
		}
	}

	matInv = matInv / fullDet;
	return std::make_tuple(matInv, true);
}


/**
 * power of a matrix, M^n
 */
template<class t_mat, class t_vec, class t_int = int>
std::tuple<t_mat, bool> pow(const t_mat& _mat, t_int pow_val)
requires is_mat<t_mat> && is_vec<t_vec> && std::is_integral_v<t_int>
{
	using t_size = decltype(_mat.size1());

	// x^0 = 1
	if(pow_val == 0)
		return std::make_tuple(unit<t_mat, t_size>(_mat.size1(), _mat.size2()), true);

	t_mat mat = _mat;

	// invert the matrix in case of negative power
	bool ok = true;
	if(pow_val < 0)
		std::tie(mat, ok) = inv<t_mat, t_vec>(mat);
	if(!ok)
		return std::make_tuple(mat, ok);

	// absolute value of the power
	const t_int pow_pos = pow_val < 0 ? -pow_val : pow_val;

	// multiply the matrix by itself |pow| times
	t_mat matpow = mat;
	for(t_int iter = 1; iter < pow_pos; ++iter)
		matpow = matpow * mat;

	return std::make_tuple(matpow, ok);
}


/**
 * gets reciprocal basis vectors |b_i> from real basis vectors |a_i> (and vice versa)
 * c: multiplicative constant (c=2*pi for physical lattices, c=1 for mathematics)
 *
 * Def.: <b_i | a_j> = c * delta(i,j)  =>
 *
 * e.g. 2d case:
 *                   ( a_1x  a_2x )
 *                   ( a_1y  a_2y )
 *
 * ( b_1x  b_1y )    (    1     0 )
 * ( b_2x  b_2y )    (    0     1 )
 *
 * B^t * A = I
 * A = B^(-t)
 */
template<class t_mat, class t_vec,
	template<class...> class t_cont_in = std::initializer_list,
	template<class...> class t_cont_out = std::vector>
t_cont_out<t_vec> recip(const t_cont_in<t_vec>& lstReal, typename t_vec::value_type c=1)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());

	const t_mat basis = create<t_mat, t_vec, t_cont_in>(lstReal);
	auto [basis_inv, bOk] = inv<t_mat, t_vec>(basis);
	basis_inv *= c;

	t_cont_out<t_vec> lstRecip;
	lstRecip.reserve(basis_inv.size1());

	for(t_size currow = 0; currow < basis_inv.size1(); ++currow)
	{
		const t_vec rowvec = row<t_mat, t_vec>(basis_inv, currow);
		lstRecip.emplace_back(std::move(rowvec));
	}

	return lstRecip;
}


/**
 * general n-dim cross product using determinant definition
 * @see https://en.wikipedia.org/wiki/Cross_product
 */
template<class t_vec, template<class...> class t_cont = std::initializer_list>
t_vec cross(const t_cont<t_vec>& vecs)
requires is_basic_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());
	using T = typename t_vec::value_type;

	// N also has to be equal to the vector size!
	const t_size N = vecs.size()+1;
	t_vec vec = zero<t_vec>(N);

	// 3-dim case
	if(N == 3 && vecs.begin()->size() == 3)
	{
		const t_vec& vec0 = *vecs.begin();
		const t_vec& vec1 = *std::next(vecs.begin(), 1);

		vec = cross<t_vec>(vec0, vec1);
	}

	// general case
	else
	{
		for(t_size iComp = 0; iComp < N; ++iComp)
		{
			std::vector<T> mat = zero<std::vector<T>>(N*N);
			mat[0*N + iComp] = T(1);

			t_size row_idx = 0;
			for(const t_vec& vec : vecs)
			{
				for(t_size col_idx = 0; col_idx < N; ++col_idx)
					mat[(row_idx+1)*N + col_idx] = vec[col_idx];
				++row_idx;
			}

			vec[iComp] = flat_det<decltype(mat)>(mat, N);
		}
	}

	return vec;
}


/**
 * intersection of plane <x|n> = d and line |org> + lam*|dir>
 * @returns [position of intersection, 0: no intersection, 1: intersection, 2: line on plane, line parameter lambda]
 * insert |x> = |org> + lam*|dir> in plane equation:
 * <org|n> + lam*<dir|n> = d
 * lam = (d - <org|n>) / <dir|n>
 *
 * @see http://mathworld.wolfram.com/Line-PlaneIntersection.html
 */
template<class t_vec>
std::tuple<t_vec, int, typename t_vec::value_type>
intersect_line_plane(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_vec& planeNorm, typename t_vec::value_type plane_d,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	// are line and plane parallel?
	const T dir_n = inner<t_vec>(lineDir, planeNorm);
	if(equals<T>(dir_n, 0, eps))
	{
		const T org_n = inner<t_vec>(lineOrg, planeNorm);
		// line on plane?
		if(equals<T>(org_n, plane_d, eps))
			return std::make_tuple(t_vec(), 2, T(0));
		// no intersection
		return std::make_tuple(t_vec(), 0, T(0));
	}

	const T org_n = inner<t_vec>(lineOrg, planeNorm);
	const T lam = (plane_d - org_n) / dir_n;

	const t_vec vecInters = lineOrg + lam*lineDir;
	return std::make_tuple(vecInters, 1, lam);
}


/**
 * intersection of plane <x|n> = d and a polygon
 * @return vertices of plane - polygon edge intersections
 */
template<class t_vec, template<class...> class t_cont>
t_cont<t_vec> intersect_plane_poly(
	const t_vec& planeNorm, typename t_vec::value_type plane_d,
	const t_cont<t_vec>& polyVerts,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
requires is_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());

	t_cont<t_vec> edgeInters;

	// intersect with each polygon edge
	for(t_size i = 0; i < polyVerts.size(); ++i)
	{
		t_size j = (i+1) % polyVerts.size();

		t_vec lineOrg = polyVerts[i];
		t_vec lineDir = polyVerts[j] - lineOrg;

		auto [pos, ty, lam] = intersect_line_plane(
			lineOrg, lineDir, planeNorm, plane_d, eps);

		if(ty == 1 && lam >= 0. && lam < 1.)
		{
			// add intersection point
			edgeInters.emplace_back(std::move(pos));
		}
		else if(ty == 2)
		{
			// line collinear to plane: add full line
			edgeInters.push_back(lineOrg);
			edgeInters.push_back(lineOrg + lineDir);
		}
	}

	return edgeInters;
}


/**
 * intersection of a sphere and a line |org> + lam*|dir>
 * @returns vector of intersections
 * insert |x> = |org> + lam*|dir> in sphere equation <x-mid | x-mid> = r^2
 *
 * @see https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection for the solution.
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_cont<t_vec>
intersect_line_sphere(
	const t_vec& lineOrg, const t_vec& _lineDir,
	const t_vec& sphereOrg, typename t_vec::value_type sphereRad,
	bool linedir_normalised = false, bool only_segment = false,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec lineDir = _lineDir;
	T lenDir = linedir_normalised ? T(1) : norm<t_vec>(lineDir);

	if(!linedir_normalised)
		lineDir /= lenDir;

	auto vecDiff = sphereOrg - lineOrg;
	auto proj = project_scalar<t_vec>(vecDiff, lineDir, true);
	auto rt = proj*proj + sphereRad*sphereRad - inner<t_vec>(vecDiff, vecDiff);

	// no intersection
	if(rt < T(0))
		return t_cont<t_vec>{};

	// one intersection
	if(equals(rt, T(0), eps))
	{
		T lam = proj/lenDir;
		if(!only_segment || (only_segment && lam >= T(0) && lam < T(1)))
			return t_cont<t_vec>{{ lineOrg + proj*lineDir }};
		return t_cont<t_vec>{};
	}

	// two intersections
	auto val = std::sqrt(rt);
	t_cont<t_vec> inters;
	inters.reserve(2);

	T lam1 = (proj + val)/lenDir;
	T lam2 = (proj - val)/lenDir;
	//std::cout << lam1 << "  " << lam2 << std::endl;

	if(!only_segment || (only_segment && lam1 >= T(0) && lam1 < T(1)))
		inters.emplace_back(lineOrg + (proj + val)*lineDir);
	if(!only_segment || (only_segment && lam2 >= T(0) && lam2 < T(1)))
		inters.emplace_back(lineOrg + (proj - val)*lineDir);

	// sort intersections by x
	std::sort(inters.begin(), inters.end(), [](const t_vec& vec1, const t_vec& vec2) -> bool
	{
		return vec1[0] < vec2[0];
	});

	return inters;
}


/**
 * intersection of two circles
 * <x-mid_{1,2} | x-mid_{1,2}> = r_{1,2}^2
 *
 * circle 1:
 * trafo to mid_1 = (0,0)
 * x^2 + y^2 = r_1^2  ->  y = +-sqrt(r_1^2 - x^2)
 *
 * circle 2:
 * (x-m_1)^2 + (y-m_2)^2 = r_2^2
 * (x-m_1)^2 + (+-sqrt(r_1^2 - x^2)-m_2)^2 = r_2^2
 *
 * @see https://mathworld.wolfram.com/Circle-CircleIntersection.html
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_cont<t_vec>
intersect_circle_circle(
	const t_vec& org1, typename t_vec::value_type r1,
	const t_vec& org2, typename t_vec::value_type r2,
	typename t_vec::value_type eps = std::sqrt(std::numeric_limits<typename t_vec::value_type>::epsilon()))
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	auto is_on_circle = [](const t_vec& org, T rad, const t_vec& pos, T eps) -> bool
	{
		T val = inner<t_vec>(org-pos, org-pos);
		return equals<T>(val, rad*rad, eps);
	};

	T m1 = org2[0] - org1[0];
	T m2 = org2[1] - org1[1];

	T r1_2 = r1*r1;
	T r2_2 = r2*r2;
	T m1_2 = m1*m1;
	T m2_2 = m2*m2;
	T m2_4 = m2_2*m2_2;

	T rt =
		+ T(2)*m2_2 * (r1_2*r2_2 + m1_2*(r1_2 + r2_2) + m2_2*(r1_2 + r2_2))
		- m2_2 * (r1_2*r1_2 + r2_2*r2_2)
		- (T(2)*m1_2*m2_4 + m1_2*m1_2*m2_2 + m2_4*m2_2);

	t_cont<t_vec> inters;
	inters.reserve(4);

	if(rt < T(0))
		return inters;

	rt = std::sqrt(rt);
	T factors = m1*(r1_2 - r2_2) + m1*m1_2 + m1*m2_2;
	T div = T(2)*(m1_2 + m2_2);

	// first intersection
	T x1 = (factors - rt) / div;
	T y1a = std::sqrt(r1_2 - x1*x1);
	T y1b = -std::sqrt(r1_2 - x1*x1);

	t_vec pos1a = m::create<t_vec>({x1, y1a}) + org1;
	t_vec pos1b = m::create<t_vec>({x1, y1b}) + org1;

	if(is_on_circle(org1, r1, pos1a, eps) && is_on_circle(org2, r2, pos1a, eps))
		inters.emplace_back(std::move(pos1a));
	if(!equals<t_vec>(pos1a, pos1b, eps) && is_on_circle(org1, r1, pos1b, eps) && is_on_circle(org2, r2, pos1b, eps))
		inters.emplace_back(std::move(pos1b));

	// second intersection
	if(!equals<T>(rt, T(0), eps))
	{
		T x2 = (factors + rt) / div;
		T y2a = std::sqrt(r1_2 - x2*x2);
		T y2b = -std::sqrt(r1_2 - x2*x2);

		t_vec pos2a = m::create<t_vec>({x2, y2a}) + org1;
		t_vec pos2b = m::create<t_vec>({x2, y2b}) + org1;

		if(is_on_circle(org1, r1, pos2a, eps) && is_on_circle(org2, r2, pos2a, eps))
			inters.emplace_back(std::move(pos2a));
		if(!equals<t_vec>(pos2a, pos2b, eps) && is_on_circle(org1, r1, pos2b, eps) && is_on_circle(org2, r2, pos2b, eps))
			inters.emplace_back(std::move(pos2b));
	}

	// sort intersections by x
	std::sort(inters.begin(), inters.end(), [](const t_vec& vec1, const t_vec& vec2) -> bool
	{
		return vec1[0] < vec2[0];
	});
	return inters;
}


/**
 * average vector or matrix
 */
template<class ty, template<class...> class t_cont = std::vector>
ty avg(const t_cont<ty>& vecs)
requires is_vec<ty> || is_mat<ty>
{
	if(vecs.size() == 0)
		return ty();

	typename ty::value_type num = 1;
	ty vec = *vecs.begin();

	auto iter = vecs.begin();
	std::advance(iter, 1);

	for(; iter != vecs.end(); std::advance(iter, 1))
	{
		vec += *iter;
		++num;
	}
	vec /= num;

	return vec;
}


/**
 * intersection of a polygon and a line
 * @returns [position of intersection, intersects?, line parameter lambda]
 */
template<class t_vec, template<class ...> class t_cont = std::vector>
std::tuple<t_vec, bool, typename t_vec::value_type>
intersect_line_poly(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_cont<t_vec>& poly)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	// middle point
	const t_vec mid = avg<t_vec, t_cont>(poly);

	// calculate polygon plane
	const t_vec vec0 = poly[0] - mid;
	const t_vec vec1 = poly[1] - mid;
	t_vec planeNorm = cross<t_vec>({vec0, vec1});
	planeNorm /= norm<t_vec>(planeNorm);
	const T planeD = inner<t_vec>(poly[0], planeNorm);

	// intersection with plane
	auto [vec, intersects, lam] = intersect_line_plane<t_vec>(lineOrg, lineDir, planeNorm, planeD);
	if(intersects != 1)
		return std::make_tuple(t_vec(), false, T(0));

	// is intersection point contained in polygon?
	const t_vec* vecFirst = &(*poly.rbegin());
	for(auto iter = poly.begin(); iter != poly.end(); std::advance(iter, 1))
	{
		const t_vec* vecSecond = &(*iter);
		const t_vec edge = *vecSecond - *vecFirst;

		// plane through edge
		t_vec edgeNorm = cross<t_vec>({edge, planeNorm});
		edgeNorm /= norm<t_vec>(edgeNorm);
		const T edgePlaneD = inner<t_vec>(*vecFirst, edgeNorm);

		// side of intersection
		const T ptEdgeD = inner<t_vec>(vec, edgeNorm);

		// outside polygon?
		if(ptEdgeD > edgePlaneD)
			return std::make_tuple(t_vec(), false, T(0));

		vecFirst = vecSecond;
	}

	// intersects with polygon
	return std::make_tuple(vec, true, lam);
}


/**
 * intersection of a polygon (transformed with a matrix) and a line
 * @returns [position of intersection, intersects?, line parameter lambda]
 */
template<class t_vec, class t_mat, template<class ...> class t_cont = std::vector>
std::tuple<t_vec, bool, typename t_vec::value_type>
intersect_line_poly(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_cont<t_vec>& _poly, const t_mat& mat)
requires is_vec<t_vec> && is_mat<t_mat>
{
	auto poly = _poly;

	// transform each vertex of the polygon
	// TODO: check for homogeneous coordinates!
	for(t_vec& vec : poly)
		vec = mat * vec;

	return intersect_line_poly<t_vec, t_cont>(lineOrg, lineDir, poly);
}


/**
 * intersection or closest points of lines |org1> + lam1*|dir1> and |org2> + lam2*|dir2>
 * @returns [nearest position 1, nearest position 2, valid?, dist, line parameter 1, line parameter 2]
 *
 * |org1> + lam1*|dir1>  =  |org2> + lam2*|dir2>
 * |org1> - |org2>  =  lam2*|dir2> - lam1*|dir1>
 * |org1> - |org2>  =  (dir2 | -dir1) * |lam2 lam1>
 * (dir2 | -dir1)^T * (|org1> - |org2>)  =  (dir2 | -dir1)^T * (dir2 | -dir1) * |lam2 lam1>
 * |lam2 lam1> = ((dir2 | -dir1)^T * (dir2 | -dir1))^(-1) * (dir2 | -dir1)^T * (|org1> - |org2>)
 *
 * @see https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
 */
template<class t_vec>
std::tuple<t_vec, t_vec, bool, typename t_vec::value_type, typename t_vec::value_type, typename t_vec::value_type>
intersect_line_line(
	const t_vec& line1Org, const t_vec& line1Dir,
	const t_vec& line2Org, const t_vec& line2Dir,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	const t_vec orgdiff = line1Org - line2Org;

	// direction matrix (symmetric)
	const T d11 = inner<t_vec>(line2Dir, line2Dir);
	const T d12 = -inner<t_vec>(line2Dir, line1Dir);
	const T d22 = inner<t_vec>(line1Dir, line1Dir);

	const T d_det = d11*d22 - d12*d12;

	// check if matrix is invertible
	if(equals<T>(d_det, 0, eps))
		return std::make_tuple(t_vec(), t_vec(), false, 0, 0, 0);

	// inverse (symmetric)
	const T d11_i = d22 / d_det;
	const T d12_i = -d12 / d_det;
	const T d22_i = d11 / d_det;

	const t_vec v1 = d11_i*line2Dir - d12_i*line1Dir;
	const t_vec v2 = d12_i*line2Dir - d22_i*line1Dir;

	const T lam2 = inner<t_vec>(v1, orgdiff);
	const T lam1 = inner<t_vec>(v2, orgdiff);

	const t_vec pos1 = line1Org + lam1*line1Dir;
	const t_vec pos2 = line2Org + lam2*line2Dir;
	const T dist = norm<t_vec>(pos2-pos1);

	return std::make_tuple(pos1, pos2, true, dist, lam1, lam2);
}


/**
 * intersection of planes <x|n1> = d1 and <x|n2> = d2
 * @returns line [org, dir, 0: no intersection, 1: intersection, 2: planes coincide]
 *
 * @see http://mathworld.wolfram.com/Plane-PlaneIntersection.html
 */
template<class t_vec>
std::tuple<t_vec, t_vec, int>
	intersect_plane_plane(
	const t_vec& plane1Norm, typename t_vec::value_type plane1_d,
	const t_vec& plane2Norm, typename t_vec::value_type plane2_d)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	//const std::size_t dim = plane1Norm.size();

	/*
	// alternate, direct calculation (TODO):
	// (n1)               ( d1 )
	// (n2) * (x,y,z)^T = ( d2 )
	//
	// N x = d  ->  R x = Q^T d  ->  back-substitute

	const t_mat N = create<t_mat, t_vec>({plane1Norm, plane2Norm}, true);
	const auto [Q, R] = qr<t_mat, t_vec>(N);

	const T d[] =
	{
		Q(0,0)*plane1_d + Q(1,0)*plane2_d,
		Q(0,1)*plane1_d + Q(1,1)*plane2_d
	};*/

	t_vec lineDir = cross<t_vec>({plane1Norm, plane2Norm});
	const T lenCross = norm<t_vec>(lineDir);

	// planes parallel or coinciding
	if(equals<T>(lenCross, 0))
	{
		const bool bCoincide = equals<T>(plane1_d, plane2_d);
		return std::make_tuple(t_vec(), t_vec(), bCoincide ? 2 : 0);
	}

	lineDir /= lenCross;

	t_vec lineOrg = - cross<t_vec>({plane1Norm, lineDir}) * plane2_d
		+ cross<t_vec>({plane2Norm, lineDir}) * plane1_d;
	lineOrg /= lenCross;

	return std::make_tuple(lineOrg, lineDir, 1);
}


/**
 * uv coordinates of a point inside a polygon defined by three vertices
 */
template<class t_vec>
t_vec poly_uv_ortho(const t_vec& vert1, const t_vec& vert2, const t_vec& vert3,
	const t_vec& uv1, const t_vec& uv2, const t_vec& uv3,
	const t_vec& _pt)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vec12 = vert2 - vert1;
	t_vec vec13 = vert3 - vert1;

	t_vec uv12 = uv2 - uv1;
	t_vec uv13 = uv3 - uv1;


	// ----------------------------------------------------
	// orthonormalisation
	const T len12 = norm<t_vec>(vec12);
	const T len13 = norm<t_vec>(vec13);
	const T lenuv12 = norm<t_vec>(uv12);
	const T lenuv13 = norm<t_vec>(uv13);
	auto vecBasis = orthonorm_sys<t_vec, std::vector, std::initializer_list>({vec12, vec13});
	auto uvBasis = orthonorm_sys<t_vec, std::vector, std::initializer_list>({uv12, uv13});
	vec12 = vecBasis[0]*len12;
	vec13 = vecBasis[1]*len13;
	uv12 = uvBasis[0]*lenuv12;
	uv13 = uvBasis[1]*lenuv13;
	// ----------------------------------------------------


	const t_vec pt = _pt - vert1;

	// project a point onto a vector and return the fraction along that vector
	auto project_lam = [](const t_vec& vec, const t_vec& vecProj) -> T
	{
		const T len = norm<t_vec>(vecProj);
		const t_vec _vecProj = vecProj / len;
		T lam = inner<t_vec>(_vecProj, vec);
		return lam / len;
	};

	T lam12 = project_lam(pt, vec12);
	T lam13 = project_lam(pt, vec13);

	// uv coordinates at specified point
	const t_vec uv_pt = uv1 + lam12*uv12 + lam13*uv13;
	return uv_pt;
}


/**
 * uv coordinates of a point inside a polygon defined by three vertices
 * (more general version than poly_uv_ortho)
 */
template<class t_mat, class t_vec>
t_vec poly_uv(const t_vec& vert1, const t_vec& vert2, const t_vec& vert3,
	const t_vec& uv1, const t_vec& uv2, const t_vec& uv3,
	const t_vec& _pt)
requires is_mat<t_mat> && is_vec<t_vec>
{
	//using T = typename t_vec::value_type;
	//using namespace m_ops;

	t_vec vec12 = vert2 - vert1;
	t_vec vec13 = vert3 - vert1;
	t_vec vecnorm = cross<t_vec>({vec12, vec13});

	// basis
	const t_mat basis = create<t_mat, t_vec>({vec12, vec13, vecnorm}, false);

	// reciprocal basis, RECI = REAL^(-T)
	const auto [basisInv, bOk] = inv<t_mat, t_vec>(basis);
	if(!bOk) return zero<t_vec>(uv1.size());

	t_vec pt = _pt - vert1;		// real pt
	pt = basisInv * pt;			// reciprocal pt

	// uv coordinates at specified point
	t_vec uv12 = uv2 - uv1;
	t_vec uv13 = uv3 - uv1;

	// pt has components in common reciprocal basis
	// assumes that both vector and uv coordinates have the same reciprocal basis
	return uv1 + pt[0]*uv12 + pt[1]*uv13;
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// 3-dim algos
// ----------------------------------------------------------------------------

/**
 * 3-dim cross product
 * @see https://en.wikipedia.org/wiki/Cross_product
 */
template<class t_vec> requires is_basic_vec<t_vec>
t_vec cross(const t_vec& vec1, const t_vec& vec2)
{
	t_vec vec;

	// only valid for 3-vectors -> use first three components
	if(vec1.size() < 3 || vec2.size() < 3)
		return vec;

	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(3);

	for(int i = 0; i < 3; ++i)
		vec[i] = vec1[(i+1)%3]*vec2[(i+2)%3] - vec1[(i+2)%3]*vec2[(i+1)%3];

	return vec;
}


/**
 * cross product matrix (3x3)
 * @see https://en.wikipedia.org/wiki/Skew-symmetric_matrix
 */
template<class t_mat, class t_vec>
t_mat skewsymmetric(const t_vec& vec)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	t_mat mat = create<t_mat>(3, 3);

	// if static matrix is larger than 3x3 (e.g. for homogeneous coordinates), initialise as identity
	if(mat.size1() > 3 || mat.size2() > 3)
		mat = unit<t_mat>(mat.size1(), mat.size2());

	mat(0,0) = 0; 		mat(0,1) = -vec[2]; 	mat(0,2) = vec[1];
	mat(1,0) = vec[2]; 	mat(1,1) = 0; 			mat(1,2) = -vec[0];
	mat(2,0) = -vec[1]; mat(2,1) = vec[0]; 		mat(2,2) = 0;

	return mat;
}


/**
 * SO(3) matrix to rotate around an axis
 * @see https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
 * @see (Arens15), p. 718 and p. 816
 * @see (Merziger06), p. 208
 */
template<class t_mat, class t_vec>
t_mat rotation(const t_vec& axis, const typename t_vec::value_type angle, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());
	using t_real = typename t_vec::value_type;

	t_real len = 1;
	if(!is_normalised)
		len = norm<t_vec>(axis);

	// ----------------------------------------------------
	// special cases: rotations around [100], [010], [001]
	if(equals(axis, create<t_vec>({len, 0, 0})))
	{
		return givens<t_mat, t_real, t_size>(axis.size(), 1, 2, angle);
		//t_real c = std::cos(angle), s = std::sin(angle);
		//return create<t_mat>({{1,0,0}, {0,c,s}, {0,-s,c}});
	}
	else if(equals(axis, create<t_vec>({0, len, 0})))
	{
		return givens<t_mat, t_real, t_size>(axis.size(), 2, 0, angle);
		//t_real c = std::cos(angle), s = std::sin(angle);
		//return create<t_mat>({{c,0,-s}, {0,1,0}, {s,0,c}});
	}
	else if(equals(axis, create<t_vec>({0, 0, len})))
	{
		return givens<t_mat, t_real, t_size>(axis.size(), 0, 1, angle);
		//t_real c = std::cos(angle), s = std::sin(angle);
		//return create<t_mat>({{c,s,0}, {-s,c,0}, {0,0,1}});
	}

	// TODO: handle 180 degrees rotation / mirroring

	// ----------------------------------------------------
	// general case
	const t_real c = std::cos(angle);
	const t_real s = std::sin(angle);

	// project along rotation axis using |v><v|
	t_mat matProj1 = projector<t_mat, t_vec>(axis/len, 1);

	// project along axis 2 in plane perpendicular to rotation axis using 1-|v><v|
	t_mat matProj2 = ortho_projector<t_mat, t_vec>(axis/len, 1) * c;

	// project along axis 3 in plane perpendicular to rotation axis and axis 2 using v_cross matrix
	t_mat matProj3 = skewsymmetric<t_mat, t_vec>(axis/len) * s;

	//std::cout << matProj1(3,3) <<  " " << matProj2(3,3) <<  " " << matProj3(3,3) << std::endl;
	// rotation in the orthogonal plane is done above by axis2*cos + axis3*sin
	t_mat matProj = matProj1 + matProj2 + matProj3;

	// if matrix is larger than 3x3 (e.g. for homogeneous coordinates), fill up with identity
	unit<t_mat>(matProj, 3,3, matProj.size1(), matProj.size2());
	return matProj;
}


/**
 * rotation axis and angle to rotate vector vec1 into vec2
 * @see O. I. Zhelezov, American Journal of Computational and Applied Mathematics 7(2), pp. 51-57 (2017), doi: 10.5923/j.ajcam.20170702.04
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_real> rotation_axis(const t_vec& vec1, const t_vec& vec2)
requires is_vec<t_vec> && is_scalar<t_real>
{
	assert(vec1.size() == vec2.size());

	// rotation axis
	t_vec axis = cross<t_vec>({ vec1, vec2 });
	const t_real lenaxis = norm<t_vec>(axis);

	// rotation angle
	const t_real& sin_angle = lenaxis;
	const t_real cos_angle = inner<t_vec>(vec1, vec2);
	const t_real angle = std::atan2(sin_angle, cos_angle);
	//std::cout << angle << " " << std::fmod(angle, pi<t_real>) << std::endl;

	axis /= lenaxis;
	return std::make_tuple(axis, angle);
}


/**
 * matrix to rotate vector vec1 into vec2
 * @see O. I. Zhelezov, American Journal of Computational and Applied Mathematics 7(2), pp. 51-57 (2017), doi: 10.5923/j.ajcam.20170702.04
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type>
t_mat rotation(const t_vec& vec1, const t_vec& vec2,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_vec<t_vec> && is_mat<t_mat> && is_scalar<t_real>
{
	assert(vec1.size() == vec2.size());
	using t_size = decltype(vec1.size());

	// rotation axis
	auto [axis, angle] = rotation_axis<t_vec, t_real>(vec1, vec2);

	// collinear vectors?
	if(equals<t_real>(angle, 0, eps))
		return unit<t_mat>(vec1.size());
	// antiparallel vectors? -> TODO
	if(equals<t_real>(std::abs(angle), pi<t_real>, eps))
	{
		t_mat mat = -unit<t_mat>(vec1.size());
		// e.g. homogeneous coordinates -> only have -1 on the first 3 diagonal elements
		for(t_size i = 3; i < std::min(mat.size1(), mat.size2()); ++i)
			mat(i,i) = 1;
		return mat;
	}

	t_mat mat = rotation<t_mat, t_vec>(axis, angle, true);
	return mat;
}


/**
 * givens rotation matrix
 * @see https://en.wikipedia.org/wiki/Givens_rotation
 */
template<class t_mat, class t_real /*= typename t_mat::value_type*/,
	typename t_size /*= decltype(t_mat{}.size1())*/>
t_mat givens(t_size N, t_size i, t_size j, t_real angle)
requires is_mat<t_mat>
{
	t_real sign = 1;

	if(j < i)
	{
		sign = -1;
		std::swap(i, j);
	}

	t_mat mat = unit<t_mat>(N, N);

	const t_real s = std::sin(angle);
	const t_real c = std::cos(angle);

	mat(i, j) = -sign*s;
	mat(j, i) = +sign*s;
	mat(i, i) = mat(j,j) = c;

	return mat;
}


/**
 * matrix to rotate n-dim vector vec1 into vec2
 * @see O. I. Zhelezov, American Journal of Computational and Applied Mathematics 7(2), pp. 51-57 (2017), doi: 10.5923/j.ajcam.20170702.04
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type>
t_mat rotation_nd(const t_vec& vec1, const t_vec& vec2)
requires is_vec<t_vec> && is_mat<t_mat>
{
	assert(vec1.size() == vec2.size());
	using t_size = decltype(vec1.size());
	const t_size dim = vec1.size();

	// rotation angle
	const t_real cos_angle = inner<t_vec>(vec1, vec2) / (norm<t_vec>(vec1) * norm<t_vec>(vec2));
	const t_real angle = std::acos(cos_angle);
	t_mat rot_01 = givens<t_mat, t_real, t_size>(dim, 0, 1, angle);

	std::vector<t_vec> vecBasis = orthonorm_sys<t_vec, std::vector>({vec1, vec2});

	t_mat matBasis, matBasis_inv;

	if(dim == 2)
	{
		matBasis = create<t_mat, t_vec>({vecBasis[0], vecBasis[1]}, true);
	}
	else if(dim == 3)
	{
		t_vec vec = cross<t_vec>({ vecBasis[0], vecBasis[1] });
		matBasis = create<t_mat, t_vec>({vecBasis[0], vecBasis[1], vec}, true);
	}
	else
	{
		// extend matBasis with linearly independent vectors until it has the full n-dim rank
		for(t_size i = 0; i < dim; ++i)
		{
			t_vec vecExt = rand<t_vec>(dim);
			vecBasis.emplace_back(std::move(vecExt));

			// TODO: check that the random vectors give full rank (i.e. det != 0)
		}

		vecBasis = orthonorm_sys<t_vec, std::vector>(vecBasis);
		vecBasis.erase(std::prev(vecBasis.end(), 2), vecBasis.end());

		matBasis = create<t_mat, t_vec>(vecBasis, true);
	}

	bool inv_ok = false;
	std::tie(matBasis_inv, inv_ok) = inv<t_mat, t_vec>(matBasis);
	//std::cout << std::boolalpha << inv_ok << std::endl;
	//matBasis_inv = trans<t_mat>(matBasis);

	return matBasis_inv * rot_01 * matBasis;
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// 3-dim algos in homogeneous coordinates
// ----------------------------------------------------------------------------

/**
 * project a homogeneous vector to screen coordinates
 * @returns [vecPersp, vecScreen]
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluProject.xml
 */
template<class t_mat, class t_vec>
std::tuple<t_vec, t_vec> hom_to_screen_coords(const t_vec& vec4,
	const t_mat& matModelView, const t_mat& matProj, const t_mat& matViewport,
	bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	// perspective trafo and divide
	t_vec vecPersp = matProj * matModelView * vec4;
	vecPersp /= vecPersp[3];

	// viewport trafo
	t_vec vec = matViewport * vecPersp;

	// flip y coordinate
	if(bFlipY) vec[1] = matViewport(1,1)*2 - vec[1];
	// flip x coordinate
	if(bFlipX) vec[0] = matViewport(0,0)*2 - vec[0];

	return std::make_tuple(vecPersp, vec);
}


/**
 * calculate world coordinates from screen coordinates
 * (vary zPlane to get the points of the z-line at constant (x,y))
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluUnProject.xml
 */
template<class t_mat, class t_vec>
t_vec hom_from_screen_coords(
	typename t_vec::value_type xScreen, typename t_vec::value_type yScreen, typename t_vec::value_type zPlane,
	const t_mat& matModelView_inv, const t_mat& matProj_inv, const t_mat& matViewport_inv,
	const t_mat* pmatViewport = nullptr, bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_vec vecScreen = create<t_vec>({xScreen, yScreen, zPlane, 1.});

	// flip y coordinate
	if(pmatViewport && bFlipY) vecScreen[1] = (*pmatViewport)(1,1)*2 - vecScreen[1];
	// flip x coordinate
	if(pmatViewport && bFlipX) vecScreen[0] = (*pmatViewport)(0,0)*2 - vecScreen[0];

	t_vec vecWorld = matModelView_inv * matProj_inv * matViewport_inv * vecScreen;

	vecWorld /= vecWorld[3];
	return vecWorld;
}


/**
 * calculate line from screen coordinates
 * @returns [pos, dir]
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluUnProject.xml
 */
template<class t_mat, class t_vec>
std::tuple<t_vec, t_vec> hom_line_from_screen_coords(
	typename t_vec::value_type xScreen, typename t_vec::value_type yScreen,
	typename t_vec::value_type z1, typename t_vec::value_type z2,
	const t_mat& matModelView_inv, const t_mat& matProj_inv, const t_mat& matViewport_inv,
	const t_mat* pmatViewport = nullptr, bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	const t_vec lineOrg = hom_from_screen_coords<t_mat, t_vec>(xScreen, yScreen, z1, matModelView_inv, matProj_inv,
		matViewport_inv, pmatViewport, bFlipY, bFlipX);
	const t_vec linePos2 = hom_from_screen_coords<t_mat, t_vec>(xScreen, yScreen, z2, matModelView_inv, matProj_inv,
		matViewport_inv, pmatViewport, bFlipY, bFlipX);

	t_vec lineDir = linePos2 - lineOrg;
	lineDir /= norm<t_vec>(lineDir);

	return std::make_tuple(lineOrg, lineDir);
}


/**
 * perspective projection matrix (homogeneous 4x4)
 * set bZ01=false for gl (near and far planes at -1 and +1), and bZ01=true for vk (planes at 0 and 1)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluPerspective.xml
 * @see https://github.com/PacktPublishing/Vulkan-Cookbook/blob/master/Library/Source%20Files/10%20Helper%20Recipes/04%20Preparing%20a%20perspective%20projection%20matrix.cpp
 * @see (Kuipers02), pp. 350-351 for a simplified version of the perspective trafo
 */
template<class t_mat>
t_mat hom_perspective(
	typename t_mat::value_type n = 0.01, typename t_mat::value_type f = 100.,
	typename t_mat::value_type fov = 0.5*pi<typename t_mat::value_type>,
	typename t_mat::value_type ratio = 3./4.,
	bool bInvZ = false, bool bZ01 = false, bool bInvY = false)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;

	const T c = 1./std::tan(0.5 * fov);
	const T n0 = bZ01 ? T(0) : n;
	const T sc = bZ01 ? T(1) : T(2);
	const T ys = bInvY ? T(-1) : T(1);
	const T zs = bInvZ ? T(-1) : T(1);

	//         ( x*c*r                           )      ( -x*c*r/z                         )
	//         ( y*c                             )      ( -y*c/z                           )
	// P * x = ( z*(n0+f)/(n-f) + w*sc*n*f/(n-f) )  =>  ( -(n0+f)/(n-f) - w/z*sc*n*f/(n-f) )
	//         ( -z                              )      ( 1                                )
	return create<t_mat>({
		c*ratio,    0.,     0.,                         0. /* t_x = 0, because it's already centred */,
		0.,         ys*c,   0.,                         0. /* t_y = 0, because it's already centred */,
		0.,         0.,     zs*(n0+f)/(n-f),            sc*n*f/(n-f),
		0.,         0.,     -zs /* persp. division */,  0.
	});
}


/**
 * simple perspective projection matrix (homogeneous 4x4), without normalising the ranges
 * @see (Kuipers02), pp. 350-351 for a simplified version of the perspective trafo
 */
template<class t_mat>
t_mat hom_perspective_no_normalisation(typename t_mat::value_type n = 0.01)
requires is_mat<t_mat>
{
	//         ( nx )    ( nx/z )
	//         ( ny )    ( ny/z )
	// P * x = ( nz ) => (  n   )
	//         (  z )    (  1   )
	return create<t_mat>({
		n,  0,  0,  0,
		0,  n,  0,  0,
		0,  0,  n,  0,
		0,  0,  1,  0
	});
}


/**
 * parallel projection matrix (homogeneous 4x4)
 * set bZ01=false for gl (near and far planes at -1 and +1), and bZ01=true for vk (planes at 0 and 1)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml
 * @see https://github.com/PacktPublishing/Vulkan-Cookbook/blob/master/Library/Source%20Files/10%20Helper%20Recipes/05%20Preparing%20an%20orthographic%20projection%20matrix.cpp
 */
template<class t_mat>
t_mat hom_parallel(
	typename t_mat::value_type n = 0.01, typename t_mat::value_type f = 100.,
	typename t_mat::value_type l = -4., typename t_mat::value_type r = 4.,
	typename t_mat::value_type b = -4., typename t_mat::value_type t = 4.,
	bool bInvZ = false, bool bZ01 = false, bool bInvY = false)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;

	const T w = r - l;
	const T h = t - b;
	const T d = n - f;

	const T sc = bZ01 ? T(1) : T(2);
	const T f0 = bZ01 ? T(0) : f;
	const T ys = bInvY ? T(-1) : T(1);
	const T zs = bInvZ ? T(-1) : T(1);

	return create<t_mat>({
		T(2)/w,   0.,         0.,       -(r+l)/w,
		0,        T(2)*ys/h,  0.,       -ys*(t+b)/h,
		0.,       0.,         sc*zs/d,   zs*(n+f0)/d,
		0.,       0.,         0.,        1.
	});
}


/**
 * parallel projection matrix (homogeneous 4x4)
 * set bZ01=false for gl (near and far planes at -1 and +1), and bZ01=true for vk (planes at 0 and 1)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml
 * @see https://github.com/PacktPublishing/Vulkan-Cookbook/blob/master/Library/Source%20Files/10%20Helper%20Recipes/05%20Preparing%20an%20orthographic%20projection%20matrix.cpp
 */
template<class t_mat>
t_mat hom_parallel_sym(
	typename t_mat::value_type n = 0.01, typename t_mat::value_type f = 100.,
	typename t_mat::value_type w = 4, typename t_mat::value_type h = 4.,
	bool bInvZ = false, bool bZ01 = false, bool bInvY = false)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;

	const T d = n - f;

	const T sc = bZ01 ? T(1) : T(2);
	const T f0 = bZ01 ? T(0) : f;
	const T ys = bInvY ? T(-1) : T(1);
	const T zs = bInvZ ? T(-1) : T(1);

	return create<t_mat>({
		T(2)/w,   0.,         0.,        0.,
		0,        T(2)*ys/h,  0.,        0.,
		0.,       0.,         sc*zs/d,   zs*(n+f0)/d,
		0.,       0.,         0.,        1.
	});
}


/**
 * viewport matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glViewport.xml
 */
template<class t_mat>
t_mat hom_viewport(typename t_mat::value_type w, typename t_mat::value_type h,
	typename t_mat::value_type n = 0, typename t_mat::value_type f = 1)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;

	return create<t_mat>({
		T(0.5)*w,  0.,        0.,            T(0.5)*w,
		0,         T(0.5)*h,  0.,            T(0.5)*h,
		0.,        0.,        T(0.5)*(f-n),  T(0.5)*(f+n),
		0.,        0.,        0.,            1.
	});
}


/**
 * translation matrix in homogeneous coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat hom_translation(t_real x, t_real y, t_real z)
requires is_mat<t_mat>
{
	return create<t_mat>({
		1.,  0.,  0.,  x,
		0.,  1.,  0.,  y,
		0.,  0.,  1.,  z,
		0.,  0.,  0.,  1.
	});
}


/**
 * translation matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_translation(const t_vec& vec)
requires is_mat<t_mat> && is_vec<t_vec>
{
	return hom_translation<t_mat, typename t_mat::value_type>(
		vec[0], vec[1], vec[2]);
}


/**
 * scaling matrix in homogeneous coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat hom_scaling(t_real x, t_real y, t_real z)
requires is_mat<t_mat>
{
	return create<t_mat>({
		x,   0.,  0.,  0.,
		0.,  y,   0.,  0.,
		0.,  0.,  z,   0.,
		0.,  0.,  0.,  1.
	});
}


/**
 * scaling matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_scaling(const t_vec& vec)
requires is_mat<t_mat> && is_vec<t_vec>
{
	return hom_scaling<t_mat, typename t_mat::value_type>(
		vec[0], vec[1], vec[2]);
}


/**
 * rotation matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_rotation(const t_vec& axis, const typename t_vec::value_type angle, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat rot = rotation<t_mat, t_vec>(axis, angle, is_normalised);

	t_mat rot_hom = unit<t_mat>(4,4);
	m::convert<t_mat, t_mat>(rot_hom, rot);

	return rot_hom;
}


/**
 * mirror matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_mirror(const t_vec& axis, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat mat = ortho_mirror_op<t_mat, t_vec>(axis, is_normalised);

	t_mat mat_hom = unit<t_mat>(4,4);
	m::convert<t_mat, t_mat>(mat_hom, mat);

	// in case the matrix is statically sized and was already larger than 4x4
	mat_hom(3, 3) = 1.;
	return mat_hom;
}


/**
 * mirror matrix in homogeneous coordinates with translation
 */
template<class t_mat, class t_vec>
t_mat hom_mirror(const t_vec& axis, const t_vec& pos, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat mirr = hom_mirror<t_mat, t_vec>(axis, is_normalised);

	t_mat offs = hom_translation<t_mat, t_vec>(pos);
	t_mat offs_inv = hom_translation<t_mat, t_vec>(-pos);
	//t_mat offs_t = trans<t_mat>(offs);

	return offs * mirr * offs_inv;
}


/**
 * "look at" matrix in homogeneous coordinates
 * @see (Sellers 2014), pp. 78-79
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluLookAt.xml
 */
template<class t_mat, class t_vec>
t_mat hom_lookat(const t_vec& pos, const t_vec& target, const t_vec& _up)
requires is_vec<t_vec> && is_mat<t_mat>
{
	//using t_real = typename t_mat::value_type;

	// right-hand system
	t_vec front = -(target - pos);
	t_vec side = cross<t_vec>({_up, front});
	t_vec up = cross<t_vec>({front, side});

	// unit vectors
	front = front / norm<t_vec>(front);
	up = up / norm<t_vec>(up);
	side = side / norm<t_vec>(side);

	// rotation matrix
	t_mat rot = unit<t_mat>(4);
	set_col<t_mat, t_vec>(rot, side, 0);
	set_col<t_mat, t_vec>(rot, up, 1);
	set_col<t_mat, t_vec>(rot, front, 2);

	// inverted rotation matrix
	t_mat rot_inv = trans<t_mat>(rot);

	/*t_vec pos_newsys = create<t_vec>({
		inner<t_vec>(pos, side),
		inner<t_vec>(pos, up),
		inner<t_vec>(pos, front)
	});*/

	// inverted translation matrix
	t_mat trans_inv = hom_translation<t_mat, t_vec>(-pos);

	// (trans * rot)^(-1)
	return rot_inv * trans_inv;
}


/**
 * shear matrix
 * @see https://en.wikipedia.org/wiki/Shear_matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type,
	typename t_size = decltype(t_mat{}.size1())>
t_mat shear(t_size ROWS, t_size COLS, t_size i, t_size j, t_real s)
requires is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(ROWS, COLS);
	mat(i,j) = s;
	return mat;
}

// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// complex algos
// ----------------------------------------------------------------------------

/**
 * SU(2) generators, pauli matrices sig_i = 2*S_i
 * @see (Arfken13), p. 110
 * @see (FUH 2021), p. 7
 */
template<class t_mat>
const t_mat& su2_matrix(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	constexpr t_cplx c0(0,0);
	constexpr t_cplx c1(1,0);
	constexpr t_cplx cI(0,1);

	static const t_mat mat[] =
	{
		create<t_mat>({{c0, c1}, { c1,  c0}}),	// x
		create<t_mat>({{c0, cI}, {-cI,  c0}}),	// y
		create<t_mat>({{c1, c0}, { c0, -c1}}),	// z
	};

	return mat[which];
}


/**
 * get a vector of pauli matrices
 * @see (Arfken13), p. 110
 */
template<class t_vec>
t_vec su2_matrices(bool bIncludeUnit = false)
requires is_basic_vec<t_vec> && is_mat<typename t_vec::value_type>
	&& is_complex<typename t_vec::value_type::value_type>
{
	using t_size = decltype(t_vec{}.size());
	using t_mat = typename t_vec::value_type;

	t_vec vec;
	vec.reserve(4);

	if(bIncludeUnit)
		vec.emplace_back(unit<t_mat>(2));
	for(t_size i = 0; i < 3; ++i)
		vec.emplace_back(su2_matrix<t_mat>(i));

	return vec;
}


/**
 * project the vector of SU(2) matrices onto a vector
 * proj = <sigma|vec>
 */
template<class t_vec, class t_mat>
t_mat proj_su2(const t_vec& vec, bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	typename t_vec::value_type len = 1;
	if(!is_normalised)
		len = norm<t_vec>(vec);

	const auto sigma = su2_matrices<std::vector<t_mat>>(false);
	return inner<std::vector<t_mat>, t_vec>(sigma, vec);
}


/**
 * SU(2) ladders
 * @see https://en.wikipedia.org/wiki/Ladder_operator
 */
template<class t_mat>
const t_mat& su2_ladder(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	constexpr t_cplx cI(0,1);
	constexpr t_cplx c05(0.5, 0);

	static const t_mat mat[] =
	{
		c05*su2_matrix<t_mat>(0) + c05*cI*su2_matrix<t_mat>(1),	// up
		c05*su2_matrix<t_mat>(0) - c05*cI*su2_matrix<t_mat>(1),	// down
	};

	return mat[which];
}


/**
 * SU(3) generators, gell-mann matrices
 * @see https://de.wikipedia.org/wiki/Gell-Mann-Matrizen
 */
template<class t_mat>
const t_mat& su3_matrix(std::size_t which)
requires is_mat<t_mat> && is_complex<typename t_mat::value_type>
{
	using t_cplx = typename t_mat::value_type;
	using t_real = typename t_cplx::value_type;
	constexpr t_cplx c0(0,0);
	constexpr t_cplx c1(1,0);
	constexpr t_cplx c2(2,0);
	constexpr t_cplx cI(0,1);
	/*constexpr*/ t_real s3 = std::sqrt(3.);

	static const t_mat mat[] =
	{
		create<t_mat>({{c0,c1,c0}, {c1,c0,c0}, {c0,c0,c0}}),			// 1
		create<t_mat>({{c0,cI,c0}, {-cI,c0,c0}, {c0,c0,c0}}),			// 2
		create<t_mat>({{c1,c0,c0}, {c0,-c1,c0}, {c0,c0,c0}}),			// 3
		create<t_mat>({{c0,c0,c1}, {c0,c0,c0}, {c1,c0,c0}}),			// 4
		create<t_mat>({{c0,c0,cI}, {c0,c0,c0}, {-cI,c0,c0}}),			// 5
		create<t_mat>({{c0,c0,c0}, {c0,c0,c1}, {c0,c1,c0}}),			// 6
		create<t_mat>({{c0,c0,c0}, {c0,c0,cI}, {c0,-cI,c0}}),			// 7
		create<t_mat>({{c1/s3,c0,c0}, {c0,c1/s3,c0}, {c0,c0,-c2/s3*c1}}),	// 8
	};

	return mat[which];
}


/**
 * crystallographic B matrix, B = 2pi * A^(-T)
 * @see https://en.wikipedia.org/wiki/Fractional_coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat B_matrix(t_real a, t_real b, t_real c, t_real _aa, t_real _bb, t_real _cc)
requires is_mat<t_mat>
{
	const t_real sc = std::sin(_cc);
	const t_real ca = std::cos(_aa);
	const t_real cb = std::cos(_bb);
	const t_real cc = std::cos(_cc);
	const t_real rr = std::sqrt(t_real{1} + t_real{2}*ca*cb*cc - (ca*ca + cb*cb + cc*cc));

	return t_real{2}*pi<t_real> * create<t_mat>({
		t_real{1}/a,			t_real{0},				t_real{0},
		t_real{-1}/a * cc/sc,		t_real{1}/b * t_real{1}/sc,		t_real{0},
		(cc*ca - cb)/(a*sc*rr), 	(cb*cc-ca)/(b*sc*rr),			sc/(c*rr)
	});
}


/**
 * crystallographic A matrix, B = 2pi * A^(-T)
 * @see https://en.wikipedia.org/wiki/Fractional_coordinates
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
t_mat A_matrix(t_real a, t_real b, t_real c, t_real _aa, t_real _bb, t_real _cc)
requires is_mat<t_mat>
{
	// |vb> derived by rotating |va>
	t_vec va = create<t_vec>({a, 0, 0});
	t_vec vb = rotation<t_mat, t_vec>(create<t_vec>({0,0,1}), _cc, 1)*va * b/a;
	t_vec vc = create<t_vec>(va.size());

	// |vc> derived using dot products:
	//   component (1) from <va|vc> = a*c*cos(_bb),
	//   component (2) from <vb|vc> = b*c*cos(_aa),
	//   component (3) from <vc|vc> = c*c
	vc[0] = std::cos(_bb)*c;
	vc[1] = (std::cos(_aa) * b*c - vb[0]*vc[0]) / vb[1];
	vc[2] = t_real{0};
	vc[2] = std::sqrt(c*c - inner<t_vec>(vc, vc));

	return create<t_mat>({
		va[0], vb[0], vc[0],
		va[1], vb[1], vc[1],
		va[2], vb[2], vc[2],
	});
}


/**
 * conjugate complex vector
 */
template<class t_vec>
t_vec conj(const t_vec& vec)
requires is_basic_vec<t_vec>
{
	using t_size = decltype(vec.size());
	const t_size N = vec.size();
	t_vec vecConj = zero<t_vec>(N);

	for(t_size iComp = 0; iComp < N; ++iComp)
	{
		if constexpr(is_complex<typename t_vec::value_type>)
			vecConj[iComp] = std::conj(vec[iComp]);
		else	// simply copy non-complex vector
			vecConj[iComp] = vec[iComp];
	}

	return vecConj;
}


/**
 * hermitian conjugate complex matrix
 */
template<class t_mat>
t_mat herm(const t_mat& mat)
requires is_basic_mat<t_mat>
{
	using t_size = decltype(mat.size1());
	t_mat mat2 = create<t_mat>(mat.size2(), mat.size1());

	for(t_size i = 0; i < mat.size1(); ++i)
		for(t_size j = 0; j < mat.size2(); ++j)
		{
			if constexpr(is_complex<typename t_mat::value_type>)
				mat2(j,i) = std::conj(mat(i,j));
			else	// simply transpose non-complex matrix
				mat2(j,i) = mat(i,j);
		}

	return mat2;
}


/**
 * bloch density operator (physically: polarisation density matrix if r = P, c = 0.5)
 *
 * eigenvector expansion of a state: |psi> = a_i |xi_i>
 * mean value of operator with mixed states:
 * <A> = p_i * <a_i|A|a_i>
 * <A> = tr( A * p_i * |a_i><a_i| )
 * <A> = tr( A * rho )
 * polarisation density matrix: rho = 0.5 * (1 + <P|sigma>)
 *
 * @see (DesktopBronstein08), Ch. 21 (Zusatzkapitel.pdf), pp. 11-12, p. 24
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9
 */
template<class t_vec, class t_mat>
t_mat bloch_density_op(const t_vec& r, typename t_vec::value_type c=0.5)
requires is_vec<t_vec> && is_mat<t_mat>
{
	return (unit<t_mat>(2,2) + proj_su2<t_vec, t_mat>(r, true)) * c;
}


/**
 * bloch vector
 * @see (DesktopBronstein08), Ch. 22 (Zusatzkapitel.pdf), p. 24
 */
template<class t_vec, class t_mat>
t_vec bloch_vector(const t_mat& state_density)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(t_vec{}.size());
	//using t_val = typename t_mat::value_type;

	const auto sigma = su2_matrices<std::vector<t_mat>>(false);
	const t_size N = sigma.size();

	t_vec bloch = create<t_vec>(N);
	for(t_size i = 0; i < N; ++i)
		bloch[i] = trace<t_mat>(state_density * sigma[i]);

	return bloch;
}

// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// quaternion algorithms
// ----------------------------------------------------------------------------

/**
 * quat1 * quat2
 * @see https://en.wikipedia.org/wiki/Quaternion#Scalar_and_vector_parts
 * @see (Kuipers02), pp. 106-110
 */
template<class t_quat>
t_quat mult(const t_quat& quat1, const t_quat& quat2)
requires m::is_quat<t_quat>
{
	/*using T = typename t_quat::value_type;
	using t_vec = std::vector<T>;

	const T r1 = quat1.real();
	//const t_vec v1 = quat1.template imag<t_vec>();
	const t_vec v1{{ quat1(1), quat1(2), quat1(3) }};
	const T r2 = quat2.real();
	//const t_vec v2 = quat2.template imag<t_vec>();
	const t_vec v2{{ quat2(1), quat2(2), quat2(3) }};

	t_quat result(r1*r2 - inner<t_vec>(v1, v2), 0, 0, 0);
	t_vec v_c = cross<t_vec>({v1, v2});

	for(int i = 0; i < 3; ++i)
		result(i+1) = v_c[i] + r1*v2[i] + r2*v1[i];

	return result;*/

	using T = typename t_quat::value_type;

	const T r1 = quat1.real();
	const T i1 = quat1.imag1();
	const T j1 = quat1.imag2();
	const T k1 = quat1.imag3();

	const T r2 = quat2.real();
	const T i2 = quat2.imag1();
	const T j2 = quat2.imag2();
	const T k2 = quat2.imag3();

	return t_quat
	{
		r1*r2 - (i1*i2 + j1*j2 + k1*k2),
		i1*r2 + r1*i2 - k1*j2 + j1*k2,
		j1*r2 + k1*i2 + r1*j2 - i1*k2,
		k1*r2 - j1*i2 + i1*j2 + r1*k2
	};
}


/**
 * conjugate quaternion
 * @see https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
 * @see also (Kuipers02), pp. 110-111
 */
template<class t_quat>
t_quat conj(const t_quat& quat)
requires m::is_quat<t_quat>
{
	return t_quat
	{
		quat.real(),
		-quat.imag1(),
		-quat.imag2(),
		-quat.imag3()
	};
}


/**
 * quat * vec
 * @see (Kuipers02), p. 127
 */
template<class t_quat, class t_vec>
t_vec mult(const t_quat& quat, const t_vec& vec)
requires m::is_quat<t_quat> && m::is_vec<t_vec>
{
	t_quat q_vec{0., vec[0], vec[1], vec[2]};

	t_quat q1 = mult<t_quat>(quat, q_vec);
	t_quat quat_rot = mult<t_quat>(q1, conj<t_quat>(quat));

	return quat_rot.template imag<t_vec>();
}


/**
 * squared quaternion norm
 * @see https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
 * @see also (Kuipers02), pp. 111-112
 */
template<class t_quat>
typename t_quat::value_type norm_sq(const t_quat& quat)
requires m::is_quat<t_quat>
{
	using t_val = typename t_quat::value_type;

	t_val r = quat.real();
	t_val i = quat.imag1();
	t_val j = quat.imag2();
	t_val k = quat.imag3();

	return r*r + i*i + j*j + k*k;
}


/**
 * quaternion norm
 * @see https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
 * @see also (Kuipers02), pp. 111-112
 */
template<class t_quat>
typename t_quat::value_type norm(const t_quat& quat)
requires m::is_quat<t_quat>
{
	return std::sqrt(norm_sq<t_quat>(quat));
}


/**
 * unit quaternion
 * @see https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
 */
template<class t_quat>
t_quat normalise(const t_quat& quat)
requires m::is_quat<t_quat>
{
	using t_val = typename t_quat::value_type;

	t_val n = norm<t_quat>(quat);
	return quat / n;
}


/**
 * inverted quaternion
 * @see https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
 * @see also (Kuipers02), p. 112
 */
template<class t_quat>
t_quat inv(const t_quat& quat)
requires is_quat<t_quat>
{
	using t_val = typename t_quat::value_type;

	t_quat quat_c = conj<t_quat>(quat);
	t_val n = norm_sq<t_quat>(quat);

	return quat_c / n;
}


/**
 * quat / quat
 * @see (DesktopBronstein08), chapter 4, equation (4.168)
 * @see (Bronstein08), chapter 4, p. 297, equation (4.115)
 * @see also (Kuipers02), p. 112
 */
template<class t_quat>
t_quat div(const t_quat& quat1, const t_quat& quat2)
requires m::is_quat<t_quat>
{
	t_quat quat2_inv = inv<t_quat>(quat2);
	return mult<t_quat>(quat1, quat2_inv);
}


/**
 * quaternion exponential
 * @see https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
 * @see (DesktopBronstein08), chapter 4, equation (4.170)
 * @see (Bronstein08), chapter 4, p. 297, equation (4.117)
 */
template<class t_quat, class t_vec>
t_quat exp(const t_quat& quat)
requires is_quat<t_quat> && is_vec<t_vec>
{
	if(equals_0<t_quat>(quat))
		return t_quat{1, 0, 0, 0};

	using t_val = typename t_quat::value_type;

	t_val r = quat.real();
	t_vec v = quat.template imag<t_vec>();
	t_val n = norm<t_vec>(v);

	t_val exp_r = std::exp(r);
	t_val ret_r = exp_r * std::cos(n);
	t_vec ret_v = exp_r * v/n*std::sin(n);

	return t_quat{ret_r, ret_v[0], ret_v[1], ret_v[2]};
}


/**
 * quaternion logarithm
 * @see https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
 */
template<class t_quat, class t_vec>
t_quat log(const t_quat& quat)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using t_val = typename t_quat::value_type;

	t_val r = quat.real();
	t_vec v = quat.template imag<t_vec>();
	t_val n_v = norm<t_vec>(v);
	t_val n_q = norm<t_quat>(quat);

	t_val ret_r = std::log(n_q);
	t_vec ret_v = v/n_v*std::acos(r/n_q);

	return t_quat{ret_r, ret_v[0], ret_v[1], ret_v[2]};
}


/**
 * quaternion power
 * @see (DesktopBronstein08), chapter 4, equation (4.180)
 * @see (Bronstein08), chapter 4, p. 298, equation (4.127)
 */
template<class t_quat, class t_vec>
t_quat pow(const t_quat& quat, typename t_quat::value_type x)
requires is_quat<t_quat> && is_vec<t_vec>
{
	t_quat l = log<t_quat, t_vec>(quat);
	return exp<t_quat, t_vec>(x * l);
}


/**
 * unit quaternion from rotation axis and angle (quaternion version of Euler formula)
 * @see https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 * @see (Bronstein08), chapter 4, pp. 301-302
 */
template<class t_quat, class t_vec, class t_real = typename t_quat::value_type>
t_quat from_rotaxis(const t_vec& vec, t_real angle)
requires is_quat<t_quat> && is_vec<t_vec>
{
	t_vec vec_norm = vec / norm<t_vec>(vec);

	t_real ret_r = std::cos(angle * t_real(0.5));
	t_vec ret_v = std::sin(angle * t_real(0.5)) * vec_norm;

	return t_quat{ret_r, ret_v[0], ret_v[1], ret_v[2]};
}


/**
 * alternate name for the above function
 */
template<class t_quat, class t_vec, class t_real = typename t_quat::value_type>
t_quat rotation(const t_vec& vec, t_real angle)
requires is_quat<t_quat> && is_vec<t_vec>
{
	return from_rotaxis<t_quat, t_vec, t_real>(vec, angle);
}


/**
 * rotation normalised axis and angle from unit quaternion (quaternion version of Euler formula)
 * @see https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 * @see (Bronstein08), chapter 4, pp. 301-302
 *
 * verifying equivalence with rodrigues formula (function rotation()) by applying a test vector x
 *
 * (a) rodrigues with (normalised) rotation axis v and angle alpha:
 * 	(|v><v| + (1-|v><v|) cos(alpha) + Px sin(alpha)) |x> =
 * 	|v><v|x> + |x> cos(alpha) - |v><v|x> cos(alpha) + v⨯x sin(alpha) =
 * 	(1 - cos(alpha)) |v><v|x> + |x> cos(alpha) + v⨯x sin(alpha)
 *
 * (b) quaternion:
 * 	(cos(alpha/2), v sin(alpha/2)) * (0, x) * (cos(alpha/2), -v sin(alpha/2)) =
 * 	(-v*x sin(alpha/2), v⨯x sin(alpha/2) + cos(alpha/2)*x) * (cos(alpha/2), -v sin(alpha/2))
 *
 * 	scalar part:
 * 		-v*x sin(alpha/2)*cos(alpha/2) + (v⨯x sin(alpha/2) + cos(alpha/2)*x)*v sin(alpha/2) =
 * 		-v*x sin(alpha/2)*cos(alpha/2) + v⨯x*v sin^2(alpha/2) + cos(alpha/2)*sin(alpha/2) x*v =
 * 		v⨯x*v sin^2(alpha/2) = 0
 *
 * 	vector part:
 * 		(v⨯x sin(alpha/2) + cos(alpha/2)*x) ⨯ (-v sin(alpha/2)) + (-v*x sin(alpha/2))*(-v sin(alpha/2)) + cos(alpha/2)*(v⨯x sin(alpha/2) + cos(alpha/2)*x) =
 * 		- v⨯x⨯v sin^2(alpha/2) - x⨯v cos(alpha/2)*sin(alpha/2) + v*x*v sin^2(alpha/2) + v⨯x sin(alpha/2)*cos(alpha/2) + x cos^2(alpha/2) =
 * 		- v⨯x⨯v sin^2(alpha/2) - 2 x⨯v cos(alpha/2)*sin(alpha/2) + v*x*v sin^2(alpha/2) + x cos^2(alpha/2) =
 * 		0.5 * (- v⨯x⨯v + v*x*v) - 0.5 *cos(alpha) * (- v⨯x⨯v + v*x*v) - x⨯v sin(alpha) + 0.5*x + 0.5*x cos(alpha) =
 * 		- 0.5*v⨯x⨯v + 0.5*v*x*v + 0.5 v⨯x⨯v cos(alpha) - 0.5*v*x*v cos(alpha) - x⨯v sin(alpha) + 0.5*x + 0.5*x cos(alpha) =
 * 		- 0.5*(x*(v*v) - v*(v*x)) + 0.5*v*x*v + 0.5 (x*(v*v) - v*(v*x)) cos(alpha) - 0.5*v*x*v cos(alpha) - x⨯v sin(alpha) + 0.5*x + 0.5*x cos(alpha) =
 * 		v*x*v - v*x*v cos(alpha) + v⨯x sin(alpha) + x cos(alpha) =
 * 		(1 - cos(alpha)) |v><v|x> + |x> cos(alpha) + v⨯x sin(alpha) ∎
 */
template<class t_quat, class t_vec, class t_real = typename t_quat::value_type>
std::tuple<t_vec, t_real> to_rotaxis(const t_quat& quat)
requires is_quat<t_quat> && is_vec<t_vec>
{
	t_real r = quat.real();
	t_vec v = quat.template imag<t_vec>();
	t_real n_q = norm<t_quat>(quat);

	t_real angle = std::acos(r/n_q);
	t_vec vec = v / (n_q * std::sin(angle));

	return std::make_tuple(vec, angle * t_real(2.));
}


/**
 * convert a quaternion to an su(2) matrix
 * @see (DesktopBronstein08), chapter 4, equations (4.163a) and (4.163b)
 * @see (Bronstein08), chapter 4, p. 296, equation (4.110a) and (4.110b)
 */
template<class t_quat, class t_vec, class t_mat_cplx,
	class t_real = typename t_quat::value_type,
	class t_cplx = typename t_mat_cplx::value_type>
t_mat_cplx to_su2(const t_quat& quat)
requires is_quat<t_quat> && is_vec<t_vec> && is_mat<t_mat_cplx>
{
	constexpr t_cplx c0(0,0);
	constexpr t_cplx c1(1,0);
	constexpr t_cplx cI(0,1);

	static const t_mat_cplx base[] =
	{
		create<t_mat_cplx>({{  c1,  c0}, {  c0,  c1 }}),	// real part
		create<t_mat_cplx>({{  c0, -cI}, { -cI,  c0 }}),	// i
		create<t_mat_cplx>({{  c0,  c1}, { -c1,  c0 }}),	// j
		create<t_mat_cplx>({{ -cI,  c0}, {  c0,  cI }}),	// k
	};

	return
		base[0] * t_cplx(quat.real()) +
		base[1] * t_cplx(quat.imag1()) +
		base[2] * t_cplx(quat.imag2()) +
		base[3] * t_cplx(quat.imag3());
}


/**
 * convert a quaternion to an so(3) matrix
 * @see (DesktopBronstein08), chapter 4, equation (4.194)
 * @see (Bronstein08), chapter 4, p. 301, equation (4.142)
 */
template<class t_quat, class t_vec, class t_mat, class t_real = typename t_quat::value_type>
t_mat to_so3(const t_quat& quat)
requires is_quat<t_quat> && is_vec<t_vec> && is_mat<t_mat>
{
	t_quat quat_conj = conj<t_quat>(quat);

	t_quat quat1 = mult<t_quat>(mult<t_quat>(quat, t_quat(0, 1, 0, 0)), quat_conj);
	t_quat quat2 = mult<t_quat>(mult<t_quat>(quat, t_quat(0, 0, 1, 0)), quat_conj);
	t_quat quat3 = mult<t_quat>(mult<t_quat>(quat, t_quat(0, 0, 0, 1)), quat_conj);

	return create<t_mat, t_vec>({
		quat1.template imag<t_vec>(),
		quat2.template imag<t_vec>(),
		quat3.template imag<t_vec>()
	});
}


/**
 * linear interpolation
 * @see (DesktopBronstein08), chapter 4, equation (4.206)
 * @see (Bronstein08), chapter 4, p. 306, equation (4.154)
 */
template<class t_quat, class t_real = typename t_quat::value_type>
t_quat lerp(const t_quat& quat1, const t_quat& quat2, t_real t)
requires is_quat<t_quat>
{
	return (t_real{1}-t)*quat1 + t*quat2;
}


/**
 * spherical linear interpolation
 * @see (DesktopBronstein08), chapter 4, equation (4.207)
 * @see (Bronstein08), chapter 4, p. 306, equation (4.155)
 */
template<class t_quat, class t_vec, class t_real = typename t_quat::value_type>
t_quat slerp(const t_quat& quat1, const t_quat& quat2, t_real t)
requires is_quat<t_quat>
{
	t_quat quat1_conj = conj<t_quat>(quat1);
	t_quat p = pow<t_quat, t_vec>(quat1_conj * quat2, t);
	return quat1 * p;
}


/**
 * quaternion to rotate vector vec1 into vec2
 */
template<class t_quat, class t_vec, class t_real = typename t_vec::value_type>
t_quat rotation(const t_vec& vec1, const t_vec& vec2, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_vec<t_vec> && is_quat<t_quat>
{
	// rotation axis
	auto [axis, angle] = rotation_axis<t_vec, t_real>(vec1, vec2);

	// collinear vectors?
	if(equals<t_real>(angle, 0, eps))
		return unit<t_quat>();
	// antiparallel vectors? -> TODO
	//if(equals<t_real>(std::abs(angle), pi<t_real>, eps))
	//	return -unit<t_quat>();

	t_quat quat = rotation<t_quat, t_vec, t_real>(axis, angle);
	return quat;
}
// ----------------------------------------------------------------------------

}
#endif
