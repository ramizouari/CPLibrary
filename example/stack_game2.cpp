#include <iostream>
#include <unordered_map>
#include <map>
#include "nt/modular_arithmetic.h"
#include <compare>
#include <tuple>
#include <random>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
#include <format>
#include <optional>
#include <future>
#include <atomic>
#include <mutex>

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
		int n = values.size(), m = O.values.size(),a=0,b=0;
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
		return (values.size()-a) <=> (O.values.size()-b);
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
		values[T.first]--;
		values[T.second]++;
		values[T.first - T.second]++;
		if (values[T.first] == 0)
			values.erase(T.first);
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
class FastState
{
	std::vector<std::pair<IntType, IntType>> values;
public:
	using AdjacentState_t = std::tuple < FastState, IntType, IntType>;
	FastState(std::vector<std::pair<IntType, IntType>> V) : values(std::move(V)) 
	{

	}
	FastState(std::vector<IntType> S)
	{
		std::sort(S.begin(), S.end());
		for (auto s : S)
			if (values.empty() || values.back().first != s)
				values.emplace_back(s, 1);
			else values.back().second++;
	}

	std::strong_ordering operator<=>(const FastState& O) const
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

	bool operator==(const FastState& O) const = default;

	std::vector<AdjacentState_t> adjacent_state() const
	{
		std::vector<AdjacentState_t> S;
		for (auto [v, _] : values) for (int k = 1; k < (v + 1) / 2; k++)
		{
			std::vector<std::pair<IntType, IntType>> adj_values;
			adj_values.reserve(values.size() + 2);
			int r;
			for (r = 0; values[r].first <= k; r++)
				adj_values.push_back(values[r]);
			if (adj_values.empty() || adj_values.back().first != k)
				adj_values.emplace_back(k, 1);
			else adj_values.back().second++;
			for (; values[r].first <= v - k; r++)
				adj_values.push_back(values[r]);
			if (adj_values.empty() || adj_values.back().first != v - k)
				adj_values.emplace_back(v - k, 1);
			else adj_values.back().second++;
			for (; r<values.size() && values[r].first <= v; r++)
				adj_values.push_back(values[r]);
			if (adj_values.back().second-- == 1)
				adj_values.pop_back();
			for (; r < values.size(); r++)
				adj_values.push_back(values[r]);
			S.emplace_back(FastState(adj_values), v, k);
		}
		return S;
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
			if(v>2)
				r += m * pow(x, v);
		return static_cast<integer>(r);
	}
};

template<typename IntType>
class std::hash<FastState<IntType>>
{
	inline static std::random_device dev;
	inline static std::mt19937_64 g = std::mt19937_64(dev());
	inline static constexpr integer M = 1e9 + 7;
	inline static std::uniform_int_distribution<integer> d = std::uniform_int_distribution<integer>(1, M - 1);
	using IK = cyclic<M>;
	std::array<IK,25> X;
public:
	hash()
	{
		for (auto& x : X)
			x = d(g);
	}
	integer operator()(const FastState<IntType>& S) const
	{
		IK r;
		int k = 0;
		for (auto [v, m] : S)
			if (v > 2)
			{
				r += v * X[k];
				r += m * X[k + 1];
				k += 2;
			}
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

template<typename IntType,typename StateType=FastState<IntType>>
class MemoizedMinimax
{
	container<StateType,result_t<IntType>> phi;

public:
	result_t<IntType> minimax(const StateType& S,bool turn)
	{
		if (phi.count(S))
		{
			const auto& [T,R, path, k] = phi.at(S);
			return { turn,T == turn ? R : -R,path,k };
		}
		IntType v=turn?-1:1;
		bool first_play = true;
		var_int r = 1;
		std::vector<std::pair<IntType,IntType>> Q;
		for (const auto& [P,m,i] : S.adjacent_state())
		{
			const auto& [_, s, path, k] = minimax(P, !turn);
			r += k;
			if (first_play || (turn && s > v || !turn && s<v))
			{
				Q.clear();
				std::copy(path.begin(), path.end(), std::back_inserter(Q));
				Q.emplace_back(m, i);
				v = s;
				first_play = false;
			}
		}
		result_t<IntType> result = result_t<IntType>{ turn, v, Q, r };
		phi[S] = result;
		return result;
	}

};

template<typename IntType,typename StateType=FastState<IntType>>
class ParallelMemoizedMinimax
{
	container<StateType, result_t<IntType>> phi;
	std::mutex lock;
	int thread = 1;
public:
	result_t<IntType> minimax(const StateType& S, bool turn,bool parallel=true)
	{
		lock.lock();
		if (phi.count(S))
		{
			const auto& [T, R, path, k] = phi.at(S);
			lock.unlock();
			return { turn,T == turn ? R : -R,path,k };
		}
		lock.unlock();
		IntType v = turn ? -1 : 1;
		bool first_play = true;
		var_int r = 1;
		std::vector<std::pair<IntType, IntType>> Q;
		if (parallel)
		{
			std::vector<std::future<result_t<IntType>>> results;
			auto Adj = S.adjacent_state();
			lock.lock();
			thread+=Adj.size();
			lock.unlock();
			for (const auto& [P, m, i] : Adj)
				results.push_back(std::async(std::launch::async, &ParallelMemoizedMinimax::minimax, this, P, !turn,false));
			
			for (auto& result : results)
				result.wait();
			lock.lock();
			thread-=Adj.size();
			lock.unlock();
			int index = 0;
			for (const auto& [P, m, i] : Adj)
			{
				const auto& [_, s, path, k] = results[index].get();
				index++;
				r += k;
				if (first_play || (turn && s > v || !turn && s < v))
				{
					Q.clear();
					std::copy(path.begin(), path.end(), std::back_inserter(Q));
					Q.emplace_back(m, i);
					v = s;
					first_play = false;
				}
			}
		}
		else for (const auto& [P, m, i] : S.adjacent_state())
		{
			const auto& [_, s, path, k] = minimax(P, !turn,false);
			r += k;
			if (first_play || (turn && s > v || !turn && s < v))
			{
				Q.clear();
				std::copy(path.begin(), path.end(), std::back_inserter(Q));
				Q.emplace_back(m, i);
				v = s;
				first_play = false;
			}
		}
		lock.lock();
		result_t<IntType> result = result_t<IntType>{ turn, v, Q, r };
		phi[S] = result;
		lock.unlock();
		return result;
	}

};



template<typename IntType>
class MemoizedPruningMinimax
{
	container<State<IntType>, result_t<IntType>> phi;
public:
	std::optional<result_t<IntType>> minimax(const State<IntType>& S, bool turn, int alpha = -1, int beta = 1)
	{
		if (phi.count(S))
		{
			const auto& [T, R, path, k] = phi.at(S);
			return result_t<IntType>{ turn,T == turn ? R : -R,path,k };
		}
		IntType v = turn ? -1 : 1;
		bool first_play = true;
		var_int r = 1;
		std::vector<std::pair<IntType, IntType>> Q;
		for (const auto& [P, m, i] : S.adjacent_state())
		{
			const auto& O = minimax(P, !turn, alpha, beta);
			if (!O.has_value())
				continue;
			const auto& [_, s, path, k] = O.value();
			r += k;
			if (alpha > beta)
				return std::nullopt;
			if (first_play || (turn && s > v || !turn && s < v))
			{
				if (turn)
					alpha = s;
				else beta = s;
				Q.clear();
				std::copy(path.begin(), path.end(), std::back_inserter(Q));
				Q.emplace_back(m, i);
				v = s;
				first_play = false;
			}
		}
		phi[S] = result_t<IntType>{ turn, v, Q, r };
		return phi[S];
	}

};

template<typename IntType>
class RandomMemoizedPruningMinimax
{
	container<State<IntType>, result_t<IntType>> phi;
public:
	std::optional<result_t<IntType>> minimax(const State<IntType>& S, bool turn, int alpha = -1, int beta = 1)
	{
		static std::random_device dev;
		static std::mt19937_64 g(dev());
		if (phi.count(S))
		{
			const auto& [T, R, path, k] = phi.at(S);
			return result_t<IntType>{ turn, T == turn ? R : -R, path, k };
		}
		IntType v = turn ? -1 : 1;
		bool first_play = true;
		var_int r = 1;
		std::vector<std::pair<IntType, IntType>> Q;
		auto R = S.adjacent_state();
		std::shuffle(R.begin(), R.end(), g);
		for (const auto& [P, m, i] : R)
		{
			const auto& O = minimax(P, !turn, alpha, beta);
			if (!O.has_value())
				continue;
			const auto& [_, s, path, k] = O.value();
			r += k;
			if (alpha > beta)
				return std::nullopt;
			if (first_play || (turn && s > v || !turn && s < v))
			{
				if (turn)
					alpha = s;
				else beta = s;
				Q.clear();
				std::copy(path.begin(), path.end(), std::back_inserter(Q));
				Q.emplace_back(m, i);
				v = s;
				first_play = false;
			}
		}
		phi[S] = result_t<IntType>{ turn, v, Q, r };
		return phi[S];
	}

};

#include <fstream>

int main()
{
	std::ofstream log("stack.txt");
	log << "n,evaluation,visited" << std::endl;
	MemoizedMinimax<char> solver;
	for (char L = 90; L >= 1; L--)
	{
		FastState state(std::vector<char>{L});
		auto [turn, v, path, k] = solver.minimax(state, true);
		log << std::format("{},{},{}", (int)L, (int)v, k.convert_to<std::string>()) << std::endl;

	}
}