//
// Created by ramizouari on 21/12/24.
//

#ifndef DIVISIBILITY_GRAPH_H
#define DIVISIBILITY_GRAPH_H
#include "graph.h"
#include "nt/number_theory.h"

namespace cp::graph {
    struct DivGraph : Graph
    {
        std::shared_ptr<abstract_factoriser> F;
        explicit DivGraph(int L,std::shared_ptr<abstract_factoriser> _f): F(std::move(_f)), Graph(L+1) {
            for (int i=1;i<=L;i++) for (auto d:F->divisors_list(i)) if (i!=d)
                connect(d,i);
            for (int i=1;i<=L;i++)
                connect(i,0);
        }

        explicit DivGraph(const std::vector<int> &E,std::shared_ptr<abstract_factoriser> _f): F(std::move(_f)), Graph(std::ranges::max(E)+1) {
            for (auto e:E) for (auto d:F->divisors_list(e)) if (e!=d)
                connect(d,e);
            for (auto e:E)
                connect(e,0);
        }
    };

    struct PrimeDivGraph : Graph
    {
        std::shared_ptr<abstract_factoriser> F;
        explicit PrimeDivGraph(int L,std::shared_ptr<abstract_factoriser> _f): F(std::move(_f)), Graph(L+1) {
            for (int i=1;i<=L;i++) for (auto p:F->prime_factors(i))
                connect(i/p,i);
        }

        explicit PrimeDivGraph(const std::vector<int>& E,std::shared_ptr<abstract_factoriser> _f): F(std::move(_f)), Graph(std::ranges::max(E)+1) {
            std::stack<int> Q;
            std::vector<bool> visited(n);
            for (auto e:E) Q.emplace(e);
            for (auto e:E) visited[e]=true;
            while (!Q.empty()) {
                auto e = Q.top();
                Q.pop();
                for (auto p:F->prime_factors(e)) {
                    connect(e/p,e);
                    if (!visited[e/p]) {
                        visited[e/p]=true;
                        Q.emplace(e/p);
                    }
                }
            }
        }
    };
}

#endif //DIVISIBILITY_GRAPH_H
