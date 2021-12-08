/**
 * hash combine test
 * @author Tobias Weber
 * @date dec-2021
 * @license see 'LICENSE' file
 */

#include <iostream>
#include <string>
#include <numeric>
#include <functional>
#include <boost/functional/hash.hpp>


template<std::size_t idx, class t_tup, class t_arr>
void assign_elem(const t_tup& tup, t_arr& hashes)
{
	using t_elem = std::tuple_element_t<idx, t_tup>;
	const t_elem& elem = std::get<idx>(tup);
	hashes[idx] = std::hash<t_elem>{}(elem);
}


template<class ...T>
std::size_t multi_hash(const T&... t)
{
	constexpr const std::size_t N = sizeof...(T);
	const std::tuple<T...> tup = std::make_tuple(t...);

	using t_hash = std::size_t;
	std::array<t_hash, N> hashes;

	auto mk_hash_arr = [&hashes, &tup]<std::size_t... idx>(const std::index_sequence<idx...>&) -> void
	{
		/*([&hashes, &tup]<std::size_t theidx>() -> void
		{
			using t_elem = std::tuple_element_t<theidx, std::tuple<T...>>;
			hashes[theidx] = std::hash<t_elem>{}(std::get<theidx>(tup));
		}<idx>(), ...);*/

		( assign_elem<idx>(tup, hashes), ...);
	};

	mk_hash_arr(std::make_index_sequence<N>());

	//for(std::size_t i=0; i<N; ++i)
	//	std::cout << std::hex << hashes[i] << std::endl;

	return std::accumulate(std::begin(hashes), std::end(hashes), t_hash(0),
		[](t_hash hash1, t_hash hash2) -> t_hash
		{
			t_hash combined_hash = hash1;
			boost::hash_combine<t_hash>(combined_hash, hash2);
			return combined_hash;
		});
}


int main()
{
	std::cout << std::hex << multi_hash<std::string>("123") << std::endl;

	std::cout << std::hex
		<< multi_hash<std::string, std::string, std::string>("123", "456", "789")
		<< std::endl;

	std::cout << std::hex << multi_hash(1, 2, 3, 4, 5) << std::endl;
	std::cout << std::hex << multi_hash(1, 2., 3.f, 4, 5) << std::endl;

	return 0;
}
