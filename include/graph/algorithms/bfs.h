//
// Created by ramizouari on 27/10/23.
//

#ifndef CPLIBRARY_BFS_H
#define CPLIBRARY_BFS_H
#include "graph/general.h"
#include <queue>

namespace cp::graph::algorithms
{
    template<typename Mapper,typename T>
    auto bfs(const AbstractGraph<T> &G,const T& u)
    {
        Mapper d;
        std::queue<T> Q;
        Q.push(u);
        d[u]=0;
        while(!Q.empty())
        {
            auto a=Q.front();
            Q.pop();
            // d[b] <= 0 means that b has not been visited yet. This is because d[b] is initialized to 0 or -inf, depending on the mapper.
            for(auto b:G.adjacentNodes(a)) if(d[b]<= 0)
            {
                d[b]=d[a]+1;
                Q.push(b);
            }
        }
        return d;
    }

    template<typename Mapper,typename T>
    auto bfs(const AbstractGraph<T> &G,const T& u,const T& v)
    {
        return bfs<Mapper,T>(G,u)[v];
    }

    template<typename Set,typename T>
    auto reachable_nodes(const AbstractGraph<T> &G,const T& u)
    {
        Set visited;
        std::queue<T> Q;
        Q.push(u);
        visited.emplace(u);
        std::vector<T> d;
        while(!Q.empty())
        {
            auto a=Q.front();
            d.push_back(a);
            Q.pop();
            // d[b] <= 0 means that b has not been visited yet. This is because d[b] is initialized to 0 or -inf, depending on the mapper.
            for(auto b:G.adjacentNodes(a)) if(!visited.count(b))
            {
                Q.push(b);
                visited.emplace(b);
            }
        }
        return d;
    }

    auto bfs(const Graph &G,int u)
    {
        std::vector<int> d(G.n);
        std::queue<int> Q;
        Q.push(u);
        d[u]=0;
        while(!Q.empty())
        {
            auto a=Q.front();
            Q.pop();
            // d[b] <= 0 means that b has not been visited yet. This is because d[b] is initialized to 0 or -inf, depending on the mapper.
            for(auto b:G.adjacentNodes(a)) if(d[b]<= 0)
                {
                    d[b]=d[a]+1;
                    Q.push(b);
                }
        }
        return d;
    }

    auto bfs(const Graph &G, int u,int v)
    {
        return bfs(G,u)[v];
    }

    auto reachable_nodes(const Graph &G,int u)
    {
        std::vector<bool> visited(G.n);
        std::queue<int> Q;
        Q.push(u);
        visited[u]=true;
        std::vector<int> d;
        while(!Q.empty())
        {
            auto a=Q.front();
            d.push_back(a);
            Q.pop();
            // d[b] <= 0 means that b has not been visited yet. This is because d[b] is initialized to 0 or -inf, depending on the mapper.
            for(auto b:G.adjacencyList[a]) if(!visited[b])
            {
                Q.push(b);
                visited[b]=true;
            }
        }
        return d;
    }

    template<typename Mapper,typename T,typename W>
    auto bfs(const AbstractWeightedGraph<T,W> &G,const T& u)
    {
        Mapper d;
        std::queue<T> Q;
        Q.push(u);
        d[u]=0;
        while(!Q.empty())
        {
            auto a=Q.front();
            Q.pop();
            // d[b] <= 0 means that b has not been visited yet. This is because d[b] is initialized to 0 or -inf, depending on the mapper.
            for(auto [b,_]:G.adjacentNodes(a)) if(d[b]<= 0)
            {
                d[b]=d[a]+1;
                Q.push(b);
            }
        }
        return d;
    }

    template<typename Mapper,typename T,typename W>
    auto bfs(const AbstractWeightedGraph<T,W> &G,const T& u,const T& v)
    {
        return bfs<Mapper,T>(G,u)[v];
    }

    template<typename Set,typename T,typename W>
    auto reachable_nodes(const AbstractWeightedGraph<T,W> &G,const T& u)
    {
        Set visited;
        std::queue<T> Q;
        Q.push(u);
        visited.emplace(u);
        std::vector<T> d;
        while(!Q.empty())
        {
            auto a=Q.front();
            d.push_back(a);
            Q.pop();
            // d[b] <= 0 means that b has not been visited yet. This is because d[b] is initialized to 0 or -inf, depending on the mapper.
            for(auto [b,_]:G.adjacentNodes(a)) if(!visited.count(b))
            {
                Q.push(b);
                visited.emplace(b);
            }
        }
        return d;
    }
}

#endif //CPLIBRARY_BFS_H
