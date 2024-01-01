//
// Created by ramizouari on 27/12/23.
//

#ifndef CPLIBRARY_FLOYD_WARSHALL_H
#define CPLIBRARY_FLOYD_WARSHALL_H
#include "algebra/order.h"
#include "graph/graph.h"

namespace cp::graph::algorithms
{
    template<typename W>
    auto floyd_warshall(const WeightedGraph<W> &G)
    {
        using H=order_closure<W>;
        std::vector<std::vector<H>> d(G.n,std::vector<H>(G.n,inf_plus));
        for(int i=0;i<G.n;i++) d[i][i]=W{};
        for(int i=0;i<G.n;i++) for(auto [j,w]:G.adjacencyList[i]) d[i][j]=std::min<H>(d[i][j],w);
        for(int k=0;k<G.n;k++) for(int i=0;i<G.n;i++) for(int j=0;j<G.n;j++)
                    d[i][j]=std::min(d[i][j],d[i][k]+d[k][j]);
        return d;
    }

    template<typename W>
    auto floyd_warshall(const WeightedGraph<W> &G,int u,int v)
    {
        return floyd_warshall(G)[u][v];
    }

    template<typename Mapper,typename T,typename W>
    auto floyd_warshall(const AbstractWeightedGraph<T,W> &G)
    {
        using H=order_closure<W>;
        Mapper d;
        for(const auto& x:G.nodes())
            d[x]=inf_plus;
        for(const auto& x:G.nodes())
            d[x][x]=W{};
        for(const auto& x:G.nodes()) for(const auto& [y,w]:G.adjacencyList[x])
                    d[x][y]=std::min(d[x][y],w);
        for(const auto& z:G.nodes()) for(const auto& x:G.nodes()) for(const auto& y:G.nodes())
                    d[x][y]=std::min(d[x][y],d[x][z]+d[z][y]);
        return d;
    }

    template<typename Mapper,typename T,typename W>
    auto floyd_warshall(const AbstractWeightedGraph<T,W> &G,const T& u,const T& v)
    {
        return floyd_warshall<Mapper,T,W>(G)[u][v];
    }
}

#endif //CPLIBRARY_FLOYD_WARSHALL_H
