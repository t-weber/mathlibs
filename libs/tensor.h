/**
 * tensor class
 * @author Tobias Weber (ident: 486fbae61c79af61ae9217361601098d7dd367380692f116995b0b81ab3b3407)
 * @date November 2021
 * @license see 'LICENSE' file
 *
 * @see references for the employed algorithms:
 * 	- (DesktopBronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its HTML version "Desktop Bronstein"].
 * 	- (Bronstein08): I. N. Bronstein et al., ISBN: 978-3-8171-2017-8 (2008) [in its paperback version].
 */

#ifndef __TENSOR_H__
#define __TENSOR_H__

#include <array>
#include <tuple>
#include <type_traits>
#include <iostream>


// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------
/**
 * multiply all function arguments
 */
template<class t_ret, class... t_args>
constexpr t_ret mult_args(const t_args&&... args) noexcept
{
	return ( args * ... );
}


template<class t_ret>
constexpr t_ret mult_args() noexcept
{
	// return 1 for rank-0 tensor (scalar)
	return t_ret{1};
}


/**
 * get the first function argument
 */
template<std::size_t i, class t_arg1>
constexpr auto&& get_arg_i(t_arg1&& arg) noexcept
{
	return arg;
}


/**
 * get the i-th function argument
 */
template<std::size_t i, class t_arg1, class... t_args>
constexpr auto&& get_arg_i(t_arg1&& arg, t_args&&... args) noexcept
{
	if constexpr(i==0)
		return arg;
	return get_arg_i<i-1>(args...);
}


/**
 * get the i-th function argument
 */
template<std::size_t i, class... t_args>
constexpr auto&& get_arg_i(t_args&&... args) noexcept
{
	return get_arg_i<i-1>(args...);
}


/**
 * concatenate two sequences
 */
template<template<std::size_t...> class t_seq = std::index_sequence,
	std::size_t ...seq1, std::size_t ...seq2>
t_seq<seq1..., seq2...>
constexpr seq_cat(const t_seq<seq1...>&, const t_seq<seq2...>&) noexcept
{
	return t_seq<seq1..., seq2...>{};
}


/**
 * first element of a sequence
 */
template<template<std::size_t...> class t_seq = std::index_sequence,
	std::size_t first, std::size_t ...seq_rest>
constexpr std::size_t seq_first(const t_seq<first, seq_rest...>&) noexcept
{
	return first;
}


/**
 * first n elements of a sequence
 */
template<template<std::size_t...> class t_seq = std::index_sequence,
	std::size_t n, std::size_t ...seq>
constexpr auto seq_first(const t_seq<seq...>&) noexcept
{
	// convert the pack into an array
	constexpr std::size_t arr[]{seq...};

	// create a new sequence with the given indices
	auto newseq = [&arr]<std::size_t ...idx>(const t_seq<idx...>&) -> auto
	{
		return t_seq<arr[idx]...>();
	};

	// take the first n indices of the array
	return newseq(std::make_index_sequence<n>());
}


/**
 * last n elements of a sequence
 */
template<template<std::size_t...> class t_seq = std::index_sequence,
	std::size_t n, std::size_t ...seq>
constexpr auto seq_last(const t_seq<seq...>&) noexcept
{
	// size of the sequence
	constexpr const std::size_t N = sizeof...(seq);

	// convert the pack into an array
	constexpr const std::size_t arr[]{seq...};

	// create a new sequence with the given indices
	auto newseq = [&arr, N]<std::size_t ...idx>(const t_seq<idx...>&) -> auto
	{
		return t_seq<arr[idx+N-n]...>();
	};

	// take the first n indices of the array
	return newseq(std::make_index_sequence<n>());
}


/**
 * remove an element from a sequence
 */
template<template<std::size_t...> class t_seq = std::index_sequence,
	std::size_t idx, std::size_t... seq>
constexpr auto seq_rm(const t_seq<seq...>& idxseq) noexcept
{
	constexpr std::size_t N = sizeof...(seq);

	auto begin = seq_first<t_seq, idx>(idxseq);
	auto end = seq_last<t_seq, N-idx-1>(idxseq);

	return seq_cat<t_seq>(begin, end);
}


/**
 * set all elements of a container to zero
 */
template<class t_cont>
void set_zero(t_cont& cont) noexcept
{
	using t_size = decltype(cont.size());
	using t_elem = std::decay_t<decltype(cont[0])>;

	for(t_size i=0; i<cont.size(); ++i)
		cont[i] = t_elem{};
}


/**
 * pick the elements with the given indices from a tuple
 */
template<std::size_t seq_offs, std::size_t ...seq>
constexpr auto pick_from_tuple(const auto& tup, const std::index_sequence<seq...>&) noexcept
{
	return std::make_tuple(std::get<seq + seq_offs>(tup)...);
}


/**
 * remove two elements from a tuple
 */
template<std::size_t idx1, std::size_t idx2, typename ...T>
constexpr auto remove_from_tuple(const std::tuple<T...>& tup) noexcept
{
	constexpr const std::size_t N = sizeof...(T);

	// before first index
	auto tup1 = pick_from_tuple<0>(tup, std::make_index_sequence<idx1>());
	// between the two indices
	auto tup2 = pick_from_tuple<idx1+1>(tup, std::make_index_sequence<idx2-idx1-1>());
	// after the second index
	auto tup3 = pick_from_tuple<idx2+1>(tup, std::make_index_sequence<N-idx2-1>());

	return std::tuple_cat(tup1, tup2, tup3);
}


/**
 * convert a container to a tuple
 */
template<class t_cont, std::size_t... seq>
auto to_tuple(const t_cont& cont, const std::index_sequence<seq...>&) noexcept
{
	return std::make_tuple(cont[seq]...);
}


/**
 * get a linear index to a multi-dimensional array with the given sizes
 */
template<class t_tup_dims, class t_tup_sizes>
constexpr std::size_t get_linear_index(const t_tup_dims& dims, const t_tup_sizes& sizes) noexcept
{
	static_assert(std::tuple_size<t_tup_dims>() == std::tuple_size<t_tup_sizes>(),
		"Wrong number of dimensions.");

	constexpr const std::size_t N = std::tuple_size<t_tup_sizes>();
	if constexpr(N == 0)
		return 0;
	else if constexpr(N == 1)
		return std::get<0>(dims);
	//else if constexpr(N == 2)
	//	return std::get<0>(dims)*std::get<1>(sizes) + std::get<1>(dims);
	else if constexpr(N > 1)
	{
		// remove first element of the sizes tuple
		auto sizes_1 = pick_from_tuple<1>(sizes, std::make_index_sequence<N-1>());
		// add a "1" to the tuple
		auto sizes_2 = std::tuple_cat(sizes_1, std::make_tuple(1));

		auto sizes_without_first = pick_from_tuple<1>(sizes_2, std::make_index_sequence<N-1>());
		auto dims_without_first = pick_from_tuple<1>(dims, std::make_index_sequence<N-1>());

		std::size_t idx =
			std::get<0>(dims) * std::get<0>(sizes_2) +
			get_linear_index(dims_without_first, sizes_without_first);

		return idx;
	}

	return 0;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
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


	/**
	 * linear element access
	 */
	constexpr t_scalar& operator[](t_size i) noexcept
	{
		return m_elems[i];
	}


	/**
	 * linear element access
	 */
	constexpr const t_scalar& operator[](t_size i) const noexcept
	{
		return m_elems[i];
	}


	/**
	 * element access
	 * TODO: make pure compile-time version without tuples
	 */
	template<class ...t_dims>
	constexpr t_scalar& operator()(const t_dims&... dims) noexcept
	{
		return operator[](::get_linear_index(
			std::forward_as_tuple(dims...),
			std::forward_as_tuple(SIZES...)));
	}


	/**
	 * element access
	 */
	template<class ...t_dims>
	constexpr const t_scalar& operator()(const t_dims&... dims) const noexcept
	{
		return operator[](::get_linear_index(
			std::forward_as_tuple(dims...),
			std::forward_as_tuple(SIZES...)));
	}


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
	constexpr t_size get_linear_index(const std::array<t_size, sizeof...(SIZES)>& arr) const
	{
		auto tup = to_tuple(arr, std::make_index_sequence<sizeof...(SIZES)>());
		return ::get_linear_index(tup, std::forward_as_tuple(SIZES...));
	}


	/**
	 * tensor contraction
	 * e.g. contract a_{ijk} over i=k: a_j = a{iji} = a{1j1} + a{2j2} + ...
	 * @see (DesktopBronstein08), ch. 4, equ. (4.75)
	 */
	template<t_size idx1, t_size idx2>
	auto contract() const noexcept
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
		res(0) = t1(0) * t2(0);
	}

	else if constexpr(rank_1 == 1 && rank_2 == 1)
	{
		for(t_size i=0; i<t1.template size<0>(); ++i)
			for(t_size j=0; j<t2.template size<0>(); ++j)
				res(i,j) = t1(i) * t2(j);
	}

	/*else if constexpr(rank_1 == 2 && rank_2 == 2)
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
// ----------------------------------------------------------------------------


#endif
