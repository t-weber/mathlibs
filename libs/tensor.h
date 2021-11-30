/**
 * statically sizes tensor class
 * @author Tobias Weber (ident: 486fbae61c79af61ae9217361601098d7dd367380692f116995b0b81ab3b3407)
 * @date November 2021
 * @license see 'LICENSE' file
 *
 * @see references for the employed algorithms:
 * 	- (DesktopBronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its HTML version "Desktop Bronstein"].
 * 	- (Bronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its paperback version].
 */

#ifndef __TENSOR_STAT_H__
#define __TENSOR_STAT_H__

#include <array>

#include "variadic_algos.h"


/**
 * tensor with static size
 */
template<class t_scalar, std::size_t ...SIZES>
class Tensor
{
public:
	using t_size = std::size_t;
	using t_cont = std::array<t_scalar, mult_args<t_size>(SIZES...)>;


public:
	constexpr Tensor(bool zero = true) noexcept
	{
		if(zero)
			set_zero<t_cont>(m_elems);
	}


	~Tensor() = default;


	/**
	 * copy constructor
	 */
	constexpr Tensor(const Tensor<t_scalar, SIZES...>& other) noexcept
	{
		operator=(other);
	}


	/**
	 * move constructor
	 */
	constexpr Tensor(Tensor<t_scalar, SIZES...>&& other) noexcept
	{
		operator=(std::forward<Tensor<t_scalar, SIZES...>&&>(other));
	}


	/**
	 * assignment
	 */
	constexpr const Tensor& operator=(const Tensor<t_scalar, SIZES...>& other) noexcept
	{
		m_elems = other.m_elems;
		return *this;
	}


	/**
	 * movement
	 */
	constexpr const Tensor& operator=(Tensor<t_scalar, SIZES...>&& other) noexcept
	{
		m_elems = std::forward<t_cont&&>(other.m_elems);
		return *this;
	}


	/**
	 * get the tensor rank
	 */
	static constexpr t_size order() noexcept
	{
		return sizeof...(SIZES);
	}


	/**
	 * get the total number of elements
	 */
	constexpr t_size size() const noexcept
	{
		return m_elems.size();
	}


	/**
	 * get the number of elements at a given position
	 */
	template<t_size i>
	static constexpr t_size size() noexcept
	{
		return get_arg_i<i>(SIZES...);
	}


	// ------------------------------------------------------------------------
	// static element access
	// ------------------------------------------------------------------------
	/**
	 * linear element access
	 */
	template<t_size I>
	constexpr t_scalar& get_lin() noexcept
	{
		return m_elems[I];
	}


	/**
	 * linear element access
	 */
	template<t_size I>
	constexpr const t_scalar& get_lin() const noexcept
	{
		return const_cast<Tensor<t_scalar, SIZES...>*>(this)->get_lin<I>();
	}


	/**
	 * element access
	 */
	template<t_size... DIMS>
	constexpr t_scalar& get() noexcept
	{
		constexpr const t_size idx =
			::get_linear_index(
				std::index_sequence<SIZES...>(),
				std::index_sequence<DIMS...>());

		return get_lin<idx>();
	}


	/**
	 * element access
	 */
	template<t_size... DIMS>
	constexpr const t_scalar& get() const noexcept
	{
		return const_cast<Tensor<t_scalar, SIZES...>*>(this)->get<DIMS...>();
	}


	/**
	 * element access
	 */
	template<t_size... DIMS>
	constexpr t_scalar& operator()() noexcept
	{
		return get<DIMS...>();
	}


	/**
	 * element access
	 */
	template<t_size... DIMS>
	constexpr const t_scalar& operator()() const noexcept
	{
		return const_cast<Tensor<t_scalar, SIZES...>*>(this)->operator()<DIMS...>();
	}
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// dynamic element access
	// ------------------------------------------------------------------------
	/**
	 * linear element access
	 */
	constexpr t_scalar& get_lin(t_size i) noexcept
	{
		return m_elems[i];
	}


	/**
	 * linear element access
	 */
	constexpr const t_scalar& get_lin(t_size i) const noexcept
	{
		return const_cast<Tensor<t_scalar, SIZES...>*>(this)->get_lin(i);
	}


	/**
	 * linear element access
	 */
	constexpr t_scalar& operator[](t_size i) noexcept
	{
		return const_cast<Tensor<t_scalar, SIZES...>*>(this)->get_lin(i);
	}


	/**
	 * linear element access
	 */
	constexpr const t_scalar& operator[](t_size i) const noexcept
	{
		return const_cast<Tensor<t_scalar, SIZES...>*>(this)->operator[](i);
	}


	/**
	 * element access
	 */
	template<class ...t_dims>
	constexpr t_scalar& get(const t_dims&... dims) noexcept
	{
		const std::array<std::size_t, sizeof...(SIZES)>
			dims_arr{{static_cast<std::size_t>(dims)...}};
		const t_size idx = ::get_linear_index<SIZES...>(dims_arr);

		return operator[](idx);

		//return operator[](::get_linear_index(
		//	std::forward_as_tuple(dims...),
		//	std::forward_as_tuple(SIZES...)));
	}


	/**
	 * element access
	 */
	template<class ...t_dims>
	constexpr const t_scalar& get(const t_dims&... dims) const noexcept
	{
		return const_cast<Tensor<t_scalar, SIZES...>*>(this)->get<t_dims...>(dims...);
	}


	/**
	 * element access
	 */
	template<class ...t_dims>
	constexpr t_scalar& operator()(const t_dims&... dims) noexcept
	{
		return get<t_dims...>(dims...);
	}


	/**
	 * element access
	 */
	template<class ...t_dims>
	constexpr const t_scalar& operator()(const t_dims&... dims) const noexcept
	{
		return const_cast<Tensor<t_scalar, SIZES...>*>(this)->operator()<t_dims...>(dims...);
	}
	// ------------------------------------------------------------------------


	/**
	 * convert a linear index to an array index
	 */
	std::array<t_size, sizeof...(SIZES)> get_array_index(t_size idx) const
	{
		constexpr const t_size order = sizeof...(SIZES);
		constexpr const t_size sizes[order]{SIZES...};

		std::array<t_size, order> arridx{};

		for(int i=(int)order-1; i>=0; --i)
		{
			arridx[i] = idx % sizes[i];
			idx /= sizes[i];
		}

		return arridx;
	}


	/**
	 * convert an array index to a linear index
	 */
	constexpr t_size get_linear_index(
		const std::array<t_size, sizeof...(SIZES)>& arr) const
	{
		return ::get_linear_index<SIZES...>(arr);

		//auto tup = to_tuple(arr, std::make_index_sequence<sizeof...(SIZES)>());
		//return ::get_linear_index(tup, std::forward_as_tuple(SIZES...));
	}


	/**
	 * tensor contraction
	 * e.g. contract a_{ijk} over i=k: a_j = a{iji} = a{1j1} + a{2j2} + ...
	 * @see (DesktopBronstein08), ch. 4, equ. (4.75)
	 */
	template<t_size idx1, t_size idx2>
	constexpr auto contract() const noexcept
	{
		// remove indices
		constexpr const t_size idx2_new = idx2 > idx1 ? idx2-1 : idx2;

		constexpr const auto seq = std::index_sequence<SIZES...>();
		constexpr const auto seq_rm1 = seq_rm<std::index_sequence, idx1>(seq);
		constexpr const auto seq_rm2 = seq_rm<std::index_sequence, idx2_new>(seq_rm1);


		// create contracted tensor
		auto t = []<t_size... idxseq>(const std::index_sequence<idxseq...>&) -> auto
		{
			return Tensor<t_scalar, idxseq...>();
		}(seq_rm2);


		// tensor types
		using t_tensor = Tensor<t_scalar, SIZES...>;
		using t_tensor_contr = std::decay_t<decltype(t)>;
		static_assert(t_tensor::size<idx1>() == t_tensor::size<idx2>(),
			"Cannot contract over indices of different dimension.");


		// calculate the contracted tensor value
		constexpr const t_size order = t_tensor::order();
		constexpr const t_size order_contr = t_tensor_contr::order();
		static_assert(order-2 == order_contr, "Wrong order of contracted tensor.");

		// create an index array to the full tensor from indices to the contracted one
		auto get_full_idx = [](const std::array<t_size, order_contr>& arr_contr)
			-> std::array<t_size, order>
		{
			std::array<t_size, order> arr{};

			for(t_size i=0; i<order_contr; ++i)
			{
				// copy the index array, leaving gaps at idx1 and idx2
				t_size idx_full = i;
				if(idx_full >= idx1) ++idx_full;
				if(idx_full >= idx2) ++idx_full;
				arr[idx_full] = arr_contr[i];
			}

			return arr;
		};

		// iterate over the components of the contracted tensor
		for(t_size idx_contr_lin=0; idx_contr_lin < t.size(); ++idx_contr_lin)
		{
			auto idx_contr_arr = t.get_array_index(idx_contr_lin);
			auto idx_full_arr = get_full_idx(idx_contr_arr);

			// iterate over the indices to contract
			for(t_size idx=0; idx<t_tensor::size<idx1>(); ++idx)
			{
				// set the two indices to contract over equal
				idx_full_arr[idx1] = idx_full_arr[idx2] = idx;

				auto idx_full_lin = get_linear_index(idx_full_arr);
				t[idx_contr_lin] += this->operator[](idx_full_lin);
			}
		}

		return t;
	}


	// ------------------------------------------------------------------------
	// operators
	/**
	 * unary +
	 */
	friend constexpr const Tensor<t_scalar, SIZES...>& operator+(
		const Tensor<t_scalar, SIZES...>& t) noexcept
	{
		return t;
	}


	/**
	 * unary -
	 */
	friend constexpr Tensor<t_scalar, SIZES...> operator-(
		const Tensor<t_scalar, SIZES...>& t) noexcept
	{
		using t_size = decltype(t.size());
		Tensor<t_scalar, SIZES...> t2;

		for(t_size i=0; i<t.size(); ++i)
			t2[i] = -t[i];

		return t2;
	}


	/**
	 * binary +
	 */
	friend constexpr Tensor<t_scalar, SIZES...> operator+(
		const Tensor<t_scalar, SIZES...>& t1, const Tensor<t_scalar, SIZES...>& t2) noexcept
	{
		using t_size = decltype(t1.size());
		Tensor<t_scalar, SIZES...> tret;

		for(t_size i=0; i<t1.size(); ++i)
			tret[i] = t1[i] + t2[i];

		return tret;
	}


	/**
	 * binary -
	 */
	friend constexpr Tensor<t_scalar, SIZES...> operator-(
		const Tensor<t_scalar, SIZES...>& t1, const Tensor<t_scalar, SIZES...>& t2) noexcept
	{
		using t_size = decltype(t1.size());
		Tensor<t_scalar, SIZES...> tret;

		for(t_size i=0; i<t1.size(); ++i)
			tret[i] = t1[i] - t2[i];

		return tret;
	}


	/**
	 * scalar multiplication
	 */
	friend constexpr Tensor<t_scalar, SIZES...> operator*(
		const Tensor<t_scalar, SIZES...>& t1, const t_scalar& s) noexcept
	{
		using t_size = decltype(t1.size());
		Tensor<t_scalar, SIZES...> tret;

		for(t_size i=0; i<t1.size(); ++i)
			tret[i] = t1[i] * s;

		return tret;
	}


	/**
	 * scalar multiplication
	 */
	friend constexpr Tensor<t_scalar, SIZES...> operator*(
		const t_scalar& s, const Tensor<t_scalar, SIZES...>& t1) noexcept
	{
		return operator*(t1, s);
	}


	/**
	 * scalar division
	 */
	friend constexpr Tensor<t_scalar, SIZES...> operator/(
		const Tensor<t_scalar, SIZES...>& t1, const t_scalar& s) noexcept
	{
		return operator*(t1, t_scalar(1)/s);
	}


	/**
	 * addition
	 */
	constexpr Tensor<t_scalar, SIZES...>& operator+=(const Tensor<t_scalar, SIZES...>& t) noexcept
	{
		for(t_size i=0; i<size(); ++i)
			operator[](i) += t[i];

		return *this;
	}


	/**
	 * subtraction
	 */
	constexpr Tensor<t_scalar, SIZES...>& operator-=(const Tensor<t_scalar, SIZES...>& t) noexcept
	{
		for(t_size i=0; i<size(); ++i)
			operator[](i) -= t[i];

		return *this;
	}


	/**
	 * scalar multiplication
	 */
	constexpr Tensor<t_scalar, SIZES...>& operator*=(const t_scalar& s) noexcept
	{
		for(t_size i=0; i<size(); ++i)
			operator[](i) *= s;

		return *this;
	}


	/**
	 * scalar division
	 */
	constexpr Tensor<t_scalar, SIZES...>& operator/=(const t_scalar& s) noexcept
	{
		return operator*=(t_scalar(1)/s);
	}
	// ------------------------------------------------------------------------


private:
	// tensor elements
	t_cont m_elems{};
};



/**
 * tensor product
 * @see (DesktopBronstein08), ch. 4, equ. (4.73a)
 */
template<class t_scalar_1, std::size_t ...SIZES_1,
	class t_scalar_2, std::size_t ...SIZES_2>
constexpr Tensor<std::common_type_t<t_scalar_1, t_scalar_2>, SIZES_1..., SIZES_2...>
tensor_prod(const Tensor<t_scalar_1, SIZES_1...>& t1, const Tensor<t_scalar_2, SIZES_2...>& t2) noexcept
{
	using t_scalar = std::common_type_t<t_scalar_1, t_scalar_2>;
	using t_tensor = Tensor<t_scalar, SIZES_1..., SIZES_2...>;
	using t_size = typename t_tensor::t_size;
	t_tensor res;

	using t_tensor_1 = std::decay_t<decltype(t1)>;
	using t_tensor_2 = std::decay_t<decltype(t2)>;
	constexpr const t_size rank_1 = t_tensor_1::order();
	constexpr const t_size rank_2 = t_tensor_2::order();

	if constexpr(rank_1 == 0 && rank_2 == 0)
	{
		res[0] = t1[0] * t2[0];
	}

	/*else if constexpr(rank_1 == 1 && rank_2 == 1)
	{
		for(t_size i=0; i<t1.template size<0>(); ++i)
			for(t_size j=0; j<t2.template size<0>(); ++j)
				res(i,j) = t1(i) * t2(j);
	}

	else if constexpr(rank_1 == 2 && rank_2 == 2)
	{
		for(t_size i=0; i<t1.template size<0>(); ++i)
			for(t_size j=0; j<t1.template size<1>(); ++j)
				for(t_size k=0; k<t2.template size<0>(); ++k)
					for(t_size l=0; l<t2.template size<1>(); ++l)
						res(i,j,k,l) = t1(i,j) * t2(k,l);
	}*/

	// general case
	else
	{
		const t_size N1 = t1.size();
		const t_size N2 = t2.size();

		for(t_size i=0; i<N1; ++i)
			for(t_size j=0; j<N2; ++j)
				res[i*N2 + j] = t1[i] * t2[j];

		// alternatively:
		//  generate rank_1 number of for loops
		//    generate rank_2 number of for loops
		//      res(rank_1 indices, rank_2 indices) = t1(rank_1 indices) * t2(rank_2 indices)
	}

	return res;
}


/**
 * tensor product
 * @see (DesktopBronstein08), ch. 4, equ. (4.73a)
 */
template<class t_scalar_1, std::size_t ...SIZES_1,
class t_scalar_2, std::size_t ...SIZES_2>
constexpr Tensor<std::common_type_t<t_scalar_1, t_scalar_2>, SIZES_1..., SIZES_2...>
operator*(const Tensor<t_scalar_1, SIZES_1...>& t1, const Tensor<t_scalar_2, SIZES_2...>& t2) noexcept
{
	return tensor_prod(t1, t2);
}


#endif
