//
// Created by ramizouari on 27/10/23.
//

#ifndef CPLIBRARY_BELLMAN_FORD_H
#define CPLIBRARY_BELLMAN_FORD_H
#include "algebra/order.h"
#include "graph/graph.h"
#include "graph/general.h"
namespace cp::graph::algorithms
{
    template<typename W>
    auto iterative_bellman_ford(const WeightedGraph<W> &G,int u,int m)
    {
        using H=order_closure<W>;
        std::vector<H> d(G.n,inf_plus);
        d[u]=W{};
        for(int i=0;i<m;i++) for(int a=0;a<G.n;a++) for(auto [b,w]:G.adjacencyList[a])
                    d[b]=std::min(d[b],d[a]+w);
        for(int i=0;i<m;i++) for(int a=0;a<G.n;a++) for(auto [b,w]:G.adjacencyList[a]) if(d[b]>d[a]+w)
                        d[b]= inf_min;
        return d;
    }

    template<typename W>
    auto bellman_ford(const WeightedGraph<W> &G,int u)
    {
        return iterative_bellman_ford(G,u,G.n-1);
    }

    template<typename W>
    auto bellman_ford(const WeightedGraph<W> &G,int u,int v)
    {
        return iterative_bellman_ford(G,u,G.n-1)[v];
    }


    template<typename Mapper,typename T,typename W>
    auto iterative_bellman_ford(const AbstractWeightedGraph<T,W> &G,const T& u,int m)
    {
        using H=order_closure<W>;
        Mapper d;
        for(const auto& x:G.nodes())
            d[x]=inf_plus;
        d[u]=W{};
        for(int i=0;i<m;i++) for(const auto& a:G.nodes()) for(const auto& [b,w]:G.adjacencyList[a])
                    d[b]=std::min(d[b],d[a]+w);
        for(int i=0;i<m;i++) for(const auto& a:G.nodes()) for(const auto& [b,w]:G.adjacencyList[a]) if(d[b]>d[a]+w)
                        d[b]= inf_min;
        return d;
    }

    template<typename Mapper,typename T,typename W>
    auto bellman_ford(const AbstractWeightedGraph<T,W> &G,const T& u)
    {
        return iterative_bellman_ford<Mapper,T,W>(G,u,G.n-1);
    }

    template<typename Mapper,typename T,typename W>
    auto bellman_ford(const AbstractWeightedGraph<T,W> &G,const T& u,const T& v)
    {
        return iterative_bellman_ford<Mapper,T,W>(G,u,G.n-1)[v];
    }
}

#endif //CPLIBRARY_BELLMAN_FORD_H
