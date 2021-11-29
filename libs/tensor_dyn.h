/**
 * dynamically sizes tensor class
 * @author Tobias Weber (ident: 486fbae61c79af61ae9217361601098d7dd367380692f116995b0b81ab3b3407)
 * @date November 2021
 * @license see 'LICENSE' file
 *
 * @see references for the employed algorithms:
 * 	- (DesktopBronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its HTML version "Desktop Bronstein"].
 * 	- (Bronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its paperback version].
 */

#ifndef __TENSOR_DYN_H__
#define __TENSOR_DYN_H__

#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>


/**
 * tensor with dynamic size
 */
template<class t_scalar, class t_size = std::size_t, template<class...> class t_cont_templ = std::vector>
class TensorDyn
{
public:
	using t_cont = t_cont_templ<t_scalar>;
	using t_cont_sizes = t_cont_templ<t_size>;


protected:
	/**
	 * get a linear index to a multi-dimensional array with the given sizes
	 */
	t_size get_linear_index(const t_cont_sizes& dims) noexcept
	{
		assert(m_sizes.size() == dims.size());

		const t_size rank = order();

		if(rank == 0)
			return 0;
		else if(rank == 1)
			return dims[0];
		else
		{
			t_size idx = 0;

			for(t_size i=0; i<rank-1; ++i)
			{
				t_size size = 1;
				for(t_size j=i+1; j<rank; ++j)
					size *= m_sizes[j];
				idx += dims[i] * size;
			}

			idx += dims[rank-1];

			return idx;
		}

		return 0;
	}


	/**
	 * convert a linear index to an array index
	 */
	t_cont_sizes get_array_index(t_size idx) const
	{
		t_cont_sizes arridx{};
		arridx.reserve(order(), 0);

		for(int i=(int)order()-1; i>=0; --i)
		{
			arridx[i] = idx % m_sizes[i];
			idx /= m_sizes[i];
		}

		return arridx;
	}


public:
	template<template<class...> class t_init = std::initializer_list>
	TensorDyn(const t_init<t_size>& sizes) noexcept
	{
		SetSizes(sizes);
	}


	TensorDyn() = default;
	~TensorDyn() = default;


	template<template<class...> class t_init = std::initializer_list>
	void SetSizes(const t_init<t_size>& sizes) noexcept
	{
		m_sizes = sizes;
		t_size m_lin_size = std::reduce(m_sizes.begin(), m_sizes.end(), 1, std::multiplies<t_size>{});
		m_elems.resize(m_lin_size, 0);
	}


	const t_cont_sizes& GetSizes() const noexcept
	{
		return m_sizes;
	}


	/**
	 * copy constructor
	 */
	TensorDyn(const TensorDyn<t_scalar, t_size, t_cont_templ>& other) noexcept
	{
		operator=(other);
	}


	/**
	 * move constructor
	 */
	TensorDyn(TensorDyn<t_scalar, t_size, t_cont_templ>&& other) noexcept
	{
		operator=(std::forward<TensorDyn<t_scalar, t_size, t_cont_templ>&&>(other));
	}


	/**
	 * assignment
	 */
	const TensorDyn& operator=(const TensorDyn<t_scalar, t_size, t_cont_templ>& other) noexcept
	{
		m_elems = other.m_elems;
		m_sizes = other.m_sizes();

		return *this;
	}


	/**
	 * movement
	 */
	const TensorDyn& operator=(TensorDyn<t_scalar, t_size, t_cont_templ>&& other) noexcept
	{
		m_elems = std::forward<t_cont&&>(other.m_elems);
		m_sizes = std::forward<t_cont_sizes&&>(other.m_sizes);

		return *this;
	}


	/**
	 * get the tensor rank
	 */
	t_size order() noexcept
	{
		return m_sizes.size();
	}


	/**
	 * get the total number of elements
	 */
	t_size size() const noexcept
	{
		return m_elems.size();
	}


	/**
	 * get the number of elements at a given position
	 */
	t_size size(t_size i) const noexcept
	{
		return m_sizes[i];
	}


	// ------------------------------------------------------------------------
	// dynamic element access
	// ------------------------------------------------------------------------
	/**
	 * linear element access
	 */
	t_scalar& get_lin(t_size i) noexcept
	{
		return m_elems[i];
	}


	/**
	 * linear element access
	 */
	const t_scalar& get_lin(t_size i) const noexcept
	{
		return const_cast<TensorDyn<t_scalar, t_size, t_cont_templ>*>(this)->get_lin(i);
	}


	/**
	 * linear element access
	 */
	t_scalar& operator[](t_size i) noexcept
	{
		return const_cast<TensorDyn<t_scalar, t_size, t_cont_templ>*>(this)->get_lin(i);
	}


	/**
	 * linear element access
	 */
	const t_scalar& operator[](t_size i) const noexcept
	{
		return const_cast<TensorDyn<t_scalar, t_size, t_cont_templ>*>(this)->operator[](i);
	}


	/**
	 * element access
	 */
	template<template<class...> class t_init = std::initializer_list>
	t_scalar& get(const t_init<t_size>& dims) noexcept
	{
		const t_size lin_idx = get_linear_index(dims);

		return operator[](lin_idx);
	}


	/**
	 * element access
	 */
	template<template<class...> class t_init = std::initializer_list>
	const t_scalar& get(const t_init<t_size>& dims) const noexcept
	{
		return const_cast<TensorDyn<t_scalar, t_size, t_cont_templ>*>(this)->get<t_init>(dims);
	}


	/**
	 * element access
	 */
	template<template<class...> class t_init = std::initializer_list>
	t_scalar& operator()(const t_init<t_size>& dims) noexcept
	{
		return get<t_init>(dims);
	}


	/**
	 * element access
	 */
	template<template<class...> class t_init = std::initializer_list>
	const t_scalar& operator()(const t_init<t_size>& dims) const noexcept
	{
		return const_cast<TensorDyn<t_scalar, t_size, t_cont_templ>*>(this)->operator()<t_init>(dims);
	}
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// operators
	/**
	 * unary +
	 */
	friend const TensorDyn<t_scalar, t_size, t_cont_templ>& operator+(
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t) noexcept
	{
		return t;
	}


	/**
	 * unary -
	 */
	friend TensorDyn<t_scalar, t_size, t_cont_templ> operator-(
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t) noexcept
	{
		TensorDyn<t_scalar, t_size, t_cont_templ> t2{t.GetSizes()};

		for(t_size i=0; i<t.size(); ++i)
			t2[i] = -t[i];

		return t2;
	}


	/**
	 * binary +
	 */
	friend TensorDyn<t_scalar, t_size, t_cont_templ> operator+(
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t1,
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t2) noexcept
	{
		assert(t1.order() == t2.order());
		TensorDyn<t_scalar, t_size, t_cont_templ> tret{t1.GetSizes()};

		for(t_size i=0; i<t1.size(); ++i)
			tret[i] = t1[i] + t2[i];

		return tret;
	}


	/**
	 * binary -
	 */
	friend TensorDyn<t_scalar, t_size, t_cont_templ> operator-(
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t1,
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t2) noexcept
	{
		assert(t1.order() == t2.order());
		TensorDyn<t_scalar, t_size, t_cont_templ> tret{t1.GetSizes()};

		for(t_size i=0; i<t1.size(); ++i)
			tret[i] = t1[i] - t2[i];

		return tret;
	}


	/**
	 * scalar multiplication
	 */
	friend TensorDyn<t_scalar, t_size, t_cont_templ> operator*(
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t1,
		const t_scalar& s) noexcept
	{
		TensorDyn<t_scalar, t_size, t_cont_templ> tret{t1.GetSizes()};

		for(t_size i=0; i<t1.size(); ++i)
			tret[i] = t1[i] * s;

		return tret;
	}


	/**
	 * scalar multiplication
	 */
	friend TensorDyn<t_scalar, t_size, t_cont_templ> operator*(
		const t_scalar& s,
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t1) noexcept
	{
		return operator*(t1, s);
	}


	/**
	 * scalar division
	 */
	friend TensorDyn<t_scalar, t_size, t_cont_templ> operator/(
		const TensorDyn<t_scalar, t_size, t_cont_templ>& t1,
		const t_scalar& s) noexcept
	{
		return operator*(t1, t_scalar(1)/s);
	}


	/**
	 * addition
	 */
	TensorDyn<t_scalar, t_size, t_cont_templ>&
	operator+=(const TensorDyn<t_scalar, t_size, t_cont_templ>& t) noexcept
	{
		assert(order() == t.order());

		for(t_size i=0; i<size(); ++i)
			operator[](i) += t[i];

		return *this;
	}


	/**
	 * subtraction
	 */
	TensorDyn<t_scalar, t_size, t_cont_templ>&
	operator-=(const TensorDyn<t_scalar, t_size, t_cont_templ>& t) noexcept
	{
		assert(order() == t.order());

		for(t_size i=0; i<size(); ++i)
			operator[](i) -= t[i];

		return *this;
	}


	/**
	 * scalar multiplication
	 */
	TensorDyn<t_scalar, t_size, t_cont_templ>&
	operator*=(const t_scalar& s) noexcept
	{
		for(t_size i=0; i<size(); ++i)
			operator[](i) *= s;

		return *this;
	}


	/**
	 * scalar division
	 */
	TensorDyn<t_scalar, t_size, t_cont_templ>&
	operator/=(const t_scalar& s) noexcept
	{
		return operator*=(t_scalar(1)/s);
	}
	// ------------------------------------------------------------------------


private:
	// tensor elements and sizes
	t_cont m_elems{};
	t_cont_sizes m_sizes{};
};


#endif