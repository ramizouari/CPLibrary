#include <iostream>
#include <unordered_map>
#include <map>
#include "modular_arithmetic.h"
#include <compare>
#include <tuple>
#include <random>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>

using var_int = boost::multiprecision::cpp_int;

template<typename K,typename V>
using container = std::map<K,V>;

constexpr std::strong_ordering reverse(const std::strong_ordering& O)
{
	return 0 <=> O;
}

template<typename IntType>
class State
{
	container<IntType,IntType> values;
public:
	using AdjacentState_t = std::tuple < State, IntType, IntType>;
	State(container<IntType, IntType> V): values(std::move(V)) {}
	State(const std::vector<IntType>& S)
	{
		for (auto s : S)
			values[s]++;
	}

	std::strong_ordering operator<=>(const State& O) const
	{
		int n = values.size(), m = O.values.size(), a = 0, b = 0;
		auto it1 = values.begin(), it2 = O.values.begin();
		while (it1 != values.end() && it1->first <= 2)
		{
			++it1;
			a++;
		}
		while (it2 != O.values.end() && it2->first <= 2)
		{
			++it2;
			b++;
		}
		for (; it1 != values.end() && it2 != O.values.end(); ++it1, ++it2)
		{
			if (it1->first != it2->first)
				return it1->first <=> it2->first;
			else if (it1->second != it2->second)
				return it2->second <=> it1->second;
		}
		return (values.size() - a) <=> (O.values.size() - b);
	}

	bool operator==(const State& O) const = default;

	std::vector<AdjacentState_t> adjacent_state() const
	{
		std::vector<AdjacentState_t> S;
		for (auto [v, _] : values) for(int k=1;k<(v+1)/2;k++)
		{
			S.emplace_back(State(values),v,k);
			std::get<0>(S.back()).values[v]--;
			if (std::get<0>(S.back()).values[v] == 0)
				std::get<0>(S.back()).values.erase(v);
			std::get<0>(S.back()).values[k]++;
			std::get<0>(S.back()).values[v - k]++;
		}
		return S;
	}

	void transform(const std::pair<IntType, IntType>& T)
	{

	}

	bool is_final_state() const
	{
		for (auto [v, m] : values)
			if (v > 2)
				return false;
		return true;
	}

	auto begin()
	{
		return values.begin();
	}

	auto begin() const
	{
		return values.begin();
	}

	auto end()
	{
		return values.end();
	}

	auto end() const
	{
		return values.end();
	}
};

template<typename IntType>
class std::hash<State<IntType>>
{
	inline static std::random_device dev;
	inline static std::mt19937_64 g = std::mt19937_64(dev());
	inline static constexpr integer M = 1e9 + 7;
	inline static std::uniform_int_distribution<integer> d = std::uniform_int_distribution<integer>(1, M - 1);
	using IK = cyclic<M>;
	IK x = d(g);
public:
	integer operator()(const State<IntType>& S) const
	{
		IK r;
		for (auto [v, m] : S)
			r += m * pow(x, v);
		return static_cast<integer>(r);
	}
};

template<typename IntType>
struct result_t
{
	using couple = std::pair<IntType, IntType>;
	bool turn;
	int evaluation;
	std::vector<couple> path;
	var_int k;
};

template<typename IntType>
class MemoizedMinimax
{
	container<State<IntType>,result_t<IntType>> phi;
public:
	result_t<IntType> minimax(const State<IntType> &S,bool turn)
	{
		if (phi.count(S))
		{
			const auto& [T,R, path, k] = phi.at(S);
			return { turn,T == turn ? R : -R,path,k };
		}
		IntType v=turn?-1:1;
		var_int r = 1;
		std::vector<std::pair<IntType,IntType>> Q;
		for (const auto& [P,m,i] : S.adjacent_state())
		{
			const auto& [_, s, path, k] = minimax(P, !turn);
			r += k;
			if (turn && s > v || !turn && s<v)
			{
				std::copy(path.begin(), path.end(), std::back_inserter(Q));
				Q.emplace_back(m, i);
				v = s;
			}
		}
		phi[S] = result_t<IntType>{ turn, v, Q, r };
		return phi[S];
	}
};


int main()
{
	MemoizedMinimax<char> solver;
	auto [turn, v, path, k] = solver.minimax(State<char>(std::vector<char>{ 60 }), true);
	std::cout << "Evaluation: " << (int)v << std::endl;
	std::reverse(path.begin(), path.end());
	std::cout << "Nodes: " << k << std::endl;
}