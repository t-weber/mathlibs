/**
 * variadic template algorithms
 * @author Tobias Weber (ident: 486fbae61c79af61ae9217361601098d7dd367380692f116995b0b81ab3b3407)
 * @date November 2021
 * @license see 'LICENSE' file
 */

#ifndef __VARIADIC_ALGOS_H__
#define __VARIADIC_ALGOS_H__

#include <array>
#include <tuple>
#include <type_traits>


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
 * get a linear index to a multi-dimensional array with the given sizes
 * (dynamically sized version)
 */
template<std::size_t ...SIZES>
constexpr std::size_t get_linear_index(
	const std::array<std::size_t, sizeof...(SIZES)>& dims) noexcept
{
	constexpr const std::size_t rank = sizeof...(SIZES);
	constexpr const std::array<std::size_t, rank> sizes{{SIZES...}};

	if constexpr(rank == 0)
		return 0;
	else if constexpr(rank == 1)
		return dims[0];
	else if constexpr(rank > 1)
	{
		std::size_t idx = 0;

		for(std::size_t i=0; i<rank-1; ++i)
		{
			std::size_t size = 1;
			for(std::size_t j=i+1; j<rank; ++j)
				size *= sizes[j];
			idx += dims[i]*size;
		}
		idx += dims[rank-1];

		return idx;
	}

	return 0;
}


/**
 * get a linear index to a multi-dimensional array with the given sizes
 * (statically sized version)
 */
template<std::size_t ...SIZES, std::size_t ...DIMS>
constexpr std::size_t get_linear_index(
	const std::index_sequence<SIZES...>&,
	const std::index_sequence<DIMS...>&) noexcept
{
	static_assert(sizeof...(SIZES) == sizeof...(DIMS),
		"Wrong number of dimensions.");

	constexpr const std::size_t rank = sizeof...(SIZES);
	constexpr const std::array<std::size_t, rank> sizes{{SIZES...}};
	constexpr const std::array<std::size_t, rank> dims{{DIMS...}};

	if constexpr(rank == 0)
		return 0;
	else if constexpr(rank == 1)
		return dims[0];
	else if constexpr(rank > 1)
	{
		std::size_t idx = 0;

		for(std::size_t i=0; i<rank-1; ++i)
		{
			std::size_t size = 1;
			for(std::size_t j=i+1; j<rank; ++j)
				size *= sizes[j];
			idx += dims[i]*size;
		}
		idx += dims[rank-1];

		return idx;
	}

	return 0;
}


// ----------------------------------------------------------------------------
// tuples
// ----------------------------------------------------------------------------
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

		// TODO: multiply sizes

		std::size_t idx =
			std::get<0>(dims) * std::get<0>(sizes_2) +
			get_linear_index(dims_without_first, sizes_without_first);

		return idx;
	}

	return 0;
}
// ----------------------------------------------------------------------------


#endif
