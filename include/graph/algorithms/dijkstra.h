//
// Created by ramizouari on 27/10/23.
//

#ifndef CPLIBRARY_DIJKSTRA_H
#define CPLIBRARY_DIJKSTRA_H
#include "graph/general.h"
#include "graph/graph.h"
#include "algebra/order.h"
#include <queue>

namespace graph::algorithms
{
    template<typename T,typename W>
    struct dijkstra_element
    {
        T node;
        order_closure<W> distance;
        dijkstra_element(T _node,order_closure<W> _distance):node(_node),distance(_distance){}
        std::strong_ordering operator<=>(const dijkstra_element& other) const
        {
            return distance<=>other.distance;
        }
    };

    template<typename W>
    auto dijkstra(const WeightedGraph<W> &G,int u)
    {
        using H=order_closure<W>;
        std::vector<H> d(G.n,inf_plus);
        d[u]=W{};
        std::priority_queue<dijkstra_element<int,W>> Q;
        Q.emplace(u,W{});
        while(!Q.empty())
        {
            auto [a,_]=Q.top();
            Q.pop();
            for(auto [b,w]:G.adjacencyList[a]) if(d[b]>d[a]+w)
                {
                    d[b]=d[a]+w;
                    Q.push({-d[b],b});
                }
        }
        return d;
    }

    template<typename W>
    auto dijkstra(const WeightedGraph<W> &G,int u,int v)
    {
        return dijkstra(G,u)[v];
    }

    template<typename Mapper,typename T,typename W>
    auto dijkstra(const AbstractWeightedGraph<T,W> &G,int u)
    {
        using H=order_closure<W>;
        Mapper d;
        for(const auto& x:G.nodes())
            d[x]=inf_plus;
        d[u]=W{};
        std::priority_queue<dijkstra_element<T,W>> Q;
        Q.emplace(u,W{});
        while(!Q.empty())
        {
            auto [a,_]=Q.top();
            Q.pop();
            for(auto [b,w]:G.adjacentNodes(a)) if(d[b]>d[a]+w)
                {
                    d[b]=d[a]+w;
                    Q.emplace(b,d[b]);
                }
        }
        return d;
    }

    template<typename Mapper,typename T,typename W>
    auto dijkstra(const AbstractWeightedGraph<T,W> &G,const T& u,const T& v)
    {
        using H=order_closure<W>;
        Mapper d;
        std::priority_queue<dijkstra_element<T,W>> Q;
        Q.emplace(u,W{});
        d[u]=W{};
        while(!Q.empty())
        {
            auto [a,q]=Q.top();
            Q.pop();
            if(d[a]<q)
                continue;
            for(auto [b,w]:G.adjacentNodes(a))
            {
                //The second condition is to deal with default values of the mapper
                if (d[b] > d[a] + w || d[b] < d[a])
                {
                    d[b] = d[a] + w;
                    Q.emplace(b, d[b]);
                }
            }
        }
        return d[v];
    }
}

#endif //CPLIBRARY_DIJKSTRA_H
