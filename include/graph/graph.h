//
// Created by ramizouari on 21/11/22.
//

#ifndef CPLIBRARY_GRAPH_H
#define CPLIBRARY_GRAPH_H

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include "union_find.h"

/**
 * @brief This the class of directed Graphs
 * @details Each Graph is a couple G=(V,E) where V={0,...,n-1} are nodes, and E is a subset of VxV
 * @param n the size of the Graph
 * */

struct Graph
{
    int n;
    std::vector<std::vector<int>> adjacencyList,reverseList;
protected:

    void topologicalSort(int r,std::vector<int> &L,std::vector<bool> &visited)
    {
        if(!visited[r])
        {
            visited[r]=true;
            for(auto s:adjacencyList[r])
                topologicalSort(s,L,visited);
            L.push_back(r);
        }
    }

    void assignComponents(int a,int u,UnionFind& C,std::vector<bool> &visited)
    {
        if(!visited[a])
        {
            visited[a]=true;
            C.connect(a,u);
            for(auto s:reverseList[a])
                assignComponents(s,u,C,visited);
        }
    }

public:
    Graph(int _n):n(_n),adjacencyList(n),reverseList(n){}
    void connect(int a,int b)
    {
        adjacencyList[a].emplace_back(b);
        reverseList[b].emplace_back(a);
    }

    std::vector<int> topologicalSort()
    {
        std::vector<bool> visited(n);
        std::vector<int> L;
        for(int i=0;i<n;i++)
            topologicalSort(i,L,visited);
        std::reverse(L.begin(),L.end());
        return L;
    }

    struct ConnectedComponentMetaData
    {
        std::vector<std::vector<int>> components;
        std::vector<int> componentId;
        UnionFind classes;
        std::vector<int> topologicalOrder;
        ConnectedComponentMetaData(std::vector<std::vector<int>> && A,std::vector<int> &&B, UnionFind&&C, std::vector<int> &&D):components(std::move(A)),
            componentId(std::move(B)), classes(std::move(C)),topologicalOrder(std::move(D)){}
    };


    ConnectedComponentMetaData getConnectedComponentsWithMetaData()
    {

        std::vector<bool> componentAssigned(n);
        UnionFind C(n);
        std::vector<int> topologicalOrder;
        for(auto l: topologicalSort())
        {
            assignComponents(l, l, C, componentAssigned);
            if(l==C.representative(l))
                topologicalOrder.push_back(l);
        }
        std::vector<bool> idAssigned(n);
        std::vector<int> componentId(n);
        std::vector<std::vector<int>> components;
        for(int i=0;i<n;i++)
        {
            auto r=C.representative(i);
            if(!idAssigned[r])
            {
                idAssigned[r]=true;
                componentId[r]=components.size();
                components.emplace_back();
            }
            componentId[i]=componentId[r];
            components[componentId[r]].push_back(i);
        }
        return ConnectedComponentMetaData(std::move(components),std::move(componentId),std::move(C),std::move(topologicalOrder));
    }

    std::pair<Graph,std::vector<int>> condensationGraph()
    {
        auto [components,componentId,classes,topologicalOrder]=getConnectedComponentsWithMetaData();
        Graph DAG(components.size());
        std::vector<std::set<int>> S(components.size());
        for(int i=0;i<n;i++) for(auto j:adjacencyList[i])
            if(!classes.equivalent(i,j) && !S[componentId[i]].contains(componentId[j]))
            {
                auto u=componentId[i],v=componentId[j];
                S[u].insert(v);
                DAG.connect(u,v);
            }
        return std::make_pair(DAG,componentId);
    }

    std::vector<std::vector<int>> getConnectedComponents()
    {
        return getConnectedComponentsWithMetaData().components;
    }
};

class PlanarGraph
{
    int n,m;
    using couple=std::pair<int,int>;
    std::vector<std::vector<std::vector<couple>>> adjacencyList,reverseList;
public:
    PlanarGraph(int n,int m):n(n),m(m),adjacencyList(n,std::vector<std::vector<couple>>(m)),reverseList(n,std::vector<std::vector<couple>>(m)){}
    void connect(int u1,int v1,int u2,int v2)
    {
        adjacencyList[u1][v1].emplace_back(u2,v2);
    }
};


/**
 * @brief This the class of directed Graphs over an ordered set
 * @details Each Graph is a couple G=(V,E) where V is an ordered set
 * @tparam OrderedSet The ordered set V
 * */
template<typename OrderedSet>
class OrderedGraph
{
    std::map<OrderedSet,std::vector<OrderedSet>> adjacencyList,reverseList;
    std::set<OrderedSet> states;

    void topologicalSort(const OrderedSet& r,std::vector<OrderedSet> &L,std::set<OrderedSet> &visited)
    {
        if(!visited.contains(r))
        {
            visited.insert(r);
            for(const auto& s:adjacencyList[r])
                topologicalSort(s,L,visited);
            L.push_back(r);
        }
    }

    void assignComponents(const OrderedSet& a,const OrderedSet& u,OrderedUnionFind<OrderedSet>& C,std::set<OrderedSet> &visited)
    {
        if(!visited.contains(a))
        {
            visited.insert(a);
            C.connect(a,u);
            for(auto s:reverseList[a])
                assignComponents(s,u,C,visited);
        }
    }

public:
    void connect(const OrderedSet&A,const OrderedSet&B)
    {
        if(!states.contains(A))
            states.insert(A);
        if(!states.contains(B))
            states.insert(B);
        adjacencyList[A].push_back(B);
        reverseList[B].push_back(A);
    }

    void addState(const OrderedSet&A)
    {
        states.insert(A);
    }

    std::vector<OrderedSet> topologicalSort()
    {
        std::set<OrderedSet> visited;
        std::vector<OrderedSet> L;
        for(auto &u:states)
            topologicalSort(u,L,visited);
        std::reverse(L.begin(),L.end());
        return L;
    }

    struct ConnectedComponentMetaData
    {
        std::vector<std::vector<OrderedSet>> components;
        std::vector<OrderedSet> componentId;
        OrderedUnionFind<OrderedSet> classes;
        std::vector<OrderedSet> topologicalOrder;
        ConnectedComponentMetaData(std::vector<std::vector<OrderedSet>> && A,std::vector<OrderedSet> &&B, OrderedUnionFind<OrderedSet>&&C,
                                   std::vector<OrderedSet>&&topologicalOrder):
        components(std::move(A)),componentId(std::move(B)), classes(std::move(C)),topologicalOrder(std::move(topologicalOrder)){}
    };

    ConnectedComponentMetaData getConnectedComponentsWithMetaData()
    {
        std::set<OrderedSet> componentAssigned;
        OrderedUnionFind<OrderedSet> C;
        std::vector<OrderedSet> topologicalOrder;
        for(auto l: topologicalSort()) {
            assignComponents(l, l, C, componentAssigned);
            if(l==C.representative(l))
                topologicalOrder.push_back(l);
        }
        std::set<OrderedSet> idAssigned;
        std::map<OrderedSet,int> componentId;
        std::vector<std::vector<OrderedSet>> components;
        for(auto &u:states)
        {
            auto r=C.representative(u);
            if(!idAssigned.contains(r))
            {
                idAssigned.insert(r);
                componentId[r]=components.size();
                components.emplace_back();
            }
            components[componentId[r]].push_back(u);
        }
        return ConnectedComponentMetaData(std::move(components),std::move(componentId),std::move(C),std::move(topologicalOrder));
    }

    std::pair<OrderedGraph,OrderedUnionFind<OrderedSet>> condensationGraph()
    {
        auto [components,componentId,classes,topologicalOrder]=getConnectedComponentsWithMetaData();
        OrderedGraph DAG(components.size());
        for(auto [u,v]:adjacencyList) if(!classes.equivalent(u,v))
            DAG.connect(classes.representative(u),classes.representative(v));
        return std::make_pair(DAG,classes);
    }

    std::vector<std::vector<OrderedSet>> getConnectedComponents()
    {
        return getConnectedComponentsWithMetaData().components;
    }
};


/**
 * @brief This the class of directed Graphs over a hashable set
 * @details Each Graph is a couple G=(V,E) where V is a hashable set
 * @tparam UnorderedSet The hashable set V
 * */
template<typename UnorderedSet>
class UnorderedGraph
{
    std::unordered_map<UnorderedSet,std::vector<UnorderedSet>> adjacencyList,reverseList;
    std::set<UnorderedSet> states;

    void topologicalSort(const UnorderedSet& r,std::vector<UnorderedSet> &L,std::unordered_set<UnorderedSet> &visited)
    {
        if(!visited.contains(r))
        {
            visited.insert(r);
            for(const auto& s:adjacencyList[r])
                topologicalSort(s,L,visited);
            L.push_back(r);
        }
    }

    void assignComponents(const UnorderedSet& a,const UnorderedSet& u,UnorderedUnionFind<UnorderedSet>& C,std::unordered_set<UnorderedSet> &visited)
    {
        if(!visited.contains(a))
        {
            visited.insert(a);
            C.connect(a,u);
            for(auto s:reverseList[a])
                assignComponents(s,u,C,visited);
        }
    }

public:
    void connect(const UnorderedSet&A,const UnorderedSet&B)
    {
        if(!states.contains(A))
            states.insert(A);
        if(!states.contains(B))
            states.insert(B);
        adjacencyList[A].push_back(B);
        reverseList[B].push_back(A);
    }

    void addState(const UnorderedSet&A)
    {
        states.insert(A);
    }

    std::vector<UnorderedSet> topologicalSort()
    {
        std::unordered_set<UnorderedSet> visited;
        std::vector<UnorderedSet> L;
        for(auto &u:states)
            topologicalSort(u,L,visited);
        std::reverse(L.begin(),L.end());
        return L;
    }

    struct ConnectedComponentMetaData
    {
        std::vector<std::vector<UnorderedSet>> components;
        std::vector<UnorderedSet> componentId;
        UnorderedUnionFind<UnorderedSet> classes;
        std::vector<UnorderedSet> topologicalOrder;
        ConnectedComponentMetaData(std::vector<std::vector<UnorderedSet>> && A,std::vector<UnorderedSet> &&B, OrderedUnionFind<UnorderedSet>&&C,
                                   std::vector<UnorderedSet> &&D):
                components(std::move(A)),componentId(std::move(B)), classes(std::move(C)),
                topologicalOrder(std::move(D)){}
    };

    ConnectedComponentMetaData getConnectedComponentsWithMetaData()
    {
        std::set<UnorderedSet> componentAssigned;
        UnorderedUnionFind<UnorderedSet> C;
        std::vector<UnorderedSet> topologicalOrder;
        for(auto l: topologicalSort()) {
            assignComponents(l, l, C, componentAssigned);
            if(l==C.representative(l))
                topologicalOrder.push_back(l);
        }
        std::unordered_set<UnorderedSet> idAssigned;
        std::unordered_map<UnorderedSet,int> componentId;
        std::vector<std::vector<UnorderedSet>> components;
        for(auto &u:states)
        {
            auto r=C.representative(u);
            if(!idAssigned.contains(r))
            {
                idAssigned.insert(r);
                componentId[r]=components.size();
                components.emplace_back();
            }
            components[componentId[r]].push_back(u);
        }
        return ConnectedComponentMetaData(std::move(components),std::move(componentId),std::move(C),std::move(topologicalOrder));
    }

    std::pair<UnorderedGraph,UnorderedUnionFind<UnorderedSet>> condensationGraph()
    {
        auto [components,componentId,classes,topologicalOrder]=getConnectedComponentsWithMetaData();
        UnorderedGraph DAG(components.size());
        for(auto [u,v]:adjacencyList) if(!classes.equivalent(u,v))
                DAG.connect(classes.representative(u),classes.representative(v));
        return std::make_pair(DAG,classes);
    }

    std::vector<std::vector<UnorderedSet>> getConnectedComponents()
    {
        return getConnectedComponentsWithMetaData().components;
    }
};


/**
 * @brief This the class of directed Graphs
 * @details Each Graph is a couple G=(V,E) where V={0,...,n-1} are nodes, and E is a subset of VxV
 * @param n the size of the Graph
 * */
template<typename Weight>
struct WeightedGraph
{
    int n;
    using AdjacentType=std::pair<int,Weight>;
    std::vector<std::vector<AdjacentType>> adjacencyList,reverseList;
protected:

    void topologicalSort(int r,std::vector<int> &L,std::vector<bool> &visited)
    {
        if(!visited[r])
        {
            visited[r]=true;
            for(auto s:adjacencyList[r])
                topologicalSort(s,L,visited);
            L.push_back(r);
        }
    }

    void assignComponents(int a,int u,UnionFind& C,std::vector<bool> &visited)
    {
        if(!visited[a])
        {
            visited[a]=true;
            C.connect(a,u);
            for(const auto& [s,_]:reverseList[a])
                assignComponents(s,u,C,visited);
        }
    }

public:
    explicit WeightedGraph(int _n):n(_n),adjacencyList(n),reverseList(n){}
    void connect(int a,int b, const Weight & w)
    {
        adjacencyList[a].emplace_back(b,w);
        reverseList[b].emplace_back(a,w);
    }

    std::vector<int> topologicalSort()
    {
        std::vector<bool> visited(n);
        std::vector<int> L;
        for(int i=0;i<n;i++)
            topologicalSort(i,L,visited);
        std::reverse(L.begin(),L.end());
        return L;
    }

    struct ConnectedComponentMetaData
    {
        std::vector<std::vector<int>> components;
        std::vector<int> componentId;
        UnionFind classes;
        std::vector<int> topologicalOrder;
        ConnectedComponentMetaData(std::vector<std::vector<int>> && A,std::vector<int> &&B, UnionFind&&C, std::vector<int> &&D):components(std::move(A)),
                                                                                                                                componentId(std::move(B)), classes(std::move(C)),topologicalOrder(std::move(D)){}
    };


    ConnectedComponentMetaData getConnectedComponentsWithMetaData()
    {

        std::vector<bool> componentAssigned(n);
        UnionFind C(n);
        std::vector<int> topologicalOrder;
        for(auto l: topologicalSort())
        {
            assignComponents(l, l, C, componentAssigned);
            if(l==C.representative(l))
                topologicalOrder.push_back(l);
        }
        std::vector<bool> idAssigned(n);
        std::vector<int> componentId(n);
        std::vector<std::vector<int>> components;
        for(int i=0;i<n;i++)
        {
            auto r=C.representative(i);
            if(!idAssigned[r])
            {
                idAssigned[r]=true;
                componentId[r]=components.size();
                components.emplace_back();
            }
            componentId[i]=componentId[r];
            components[componentId[r]].push_back(i);
        }
        return ConnectedComponentMetaData(std::move(components),std::move(componentId),std::move(C),std::move(topologicalOrder));
    }

    std::pair<Graph,std::vector<int>> condensationGraph()
    {
        auto [components,componentId,classes,topologicalOrder]=getConnectedComponentsWithMetaData();
        Graph DAG(components.size());
        std::vector<std::set<int>> S(components.size());
        for(int i=0;i<n;i++) for(auto j:adjacencyList[i])
                if(!classes.equivalent(i,j) && !S[componentId[i]].contains(componentId[j]))
                {
                    auto u=componentId[i],v=componentId[j];
                    S[u].insert(v);
                    DAG.connect(u,v);
                }
        return std::make_pair(DAG,componentId);
    }

    std::vector<std::vector<int>> getConnectedComponents()
    {
        return getConnectedComponentsWithMetaData().components;
    }
};


#endif //CPLIBRARY_GRAPH_H
