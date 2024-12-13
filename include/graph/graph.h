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
#include <span>
#include "union_find.h"
#include "general.h"

/**
 * @brief This the class of directed Graphs
 * @details Each Graph is a couple G=(V,E) where V={0,...,n-1} are nodes, and E is a subset of VxV
 * @param n the size of the Graph
 * */

namespace cp::graph
{
    struct Graph : AbstractGraph<int>
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

            // Decomposition of strongly connected components
            struct SCDec {
                int components; // Number of components
                std::vector<int> componentId,// A map that gives the id of the strongly connected component of a node
                                 topologicalOrder; // A topological order of the components
            };

            void assign(int u, std::vector<bool> &assigned,std::vector<int> & cmpId,int id)
            {
                if (assigned[u]) return;
                cmpId[u]=id;
                assigned[u] = true;
                for (int v : reverseList[u]) assign(v,assigned,cmpId,id);
            }

            SCDec getSCComponents()
            {
                std::vector<bool> hasId(n);
                std::vector<int> topoOrder = topologicalSort(),cmpId(n);
                int components=0;
                for (auto i :topoOrder) if (!hasId[i])
                    assign(i,hasId,cmpId,components++);
                return {components,cmpId, topoOrder};
            }

            std::pair<Graph,std::vector<int>> condensationGraph()
            {
                auto [components,componentId,topologicalOrder]=getSCComponents();
                Graph DAG(components);
                std::vector<std::set<int>> S(components);
                for(int i=0;i<n;i++) for(auto j:adjacencyList[i]) if(componentId[i]!=componentId[j] && !S[componentId[i]].contains(componentId[j]))
                {
                    auto u=componentId[i],v=componentId[j];
                    S[u].insert(v);
                    DAG.connect(u,v);
                }
                return std::make_pair(DAG,componentId);
            }


        int size() const override
        {
            return n;
        }

        view_or_value<int> adjacentNodes(const int&u, bool direction) const override
        {
            if(direction)
                return std::span<const int>(adjacencyList[u].data(),adjacencyList[u].size());
            return std::span<const int>(reverseList[u].data(),reverseList[u].size());
        }

        view_or_value<int> adjacentNodes(const int&u) const override
        {
            return adjacentNodes(u,true);
        }

        view_or_value<int> nodes() const override
        {
            std::vector<int> A(n);
            for(int i=0;i<n;i++)
                A[i]=i;
            return A;
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

        // Decomposition of strongly connected components
        struct SCDec {
            int components; // Number of components
            std::map<OrderedSet,int> componentId;// A map that gives the id of the strongly connected component of a node
            std::vector<OrderedSet> topologicalOrder; // A topological order of the components
        };

        void assign(const OrderedSet &u, std::set<OrderedSet> &assigned,std::map<OrderedSet,int> & cmpId,int id)
        {
            if (assigned.contains(u)) return;
            cmpId[u]=id;
            assigned.insert(u);
            for (const auto& v : reverseList[u]) assign(v,assigned,cmpId,id);
        }

        SCDec getSCComponents()
        {
            std::set<OrderedSet> hasId;
            std::vector<OrderedSet> topoOrder = topologicalSort();
            std::map<OrderedSet,int> cmpId;
            int components=0;
            for (const auto& u :topoOrder) if (!hasId.contains(u))
                assign(u,hasId,cmpId,components++);
            return {components,cmpId, topoOrder};
        }

        std::pair<Graph,std::map<OrderedSet,int>> condensationGraph()
        {
            auto [components,componentId,topologicalOrder]=getSCComponents();
            Graph DAG(components);
            std::vector<std::set<int>> S(components);
            for(int i=0;i<states.size();i++) for(auto j:adjacencyList[i]) if(componentId[i]!=componentId[j] && !S[componentId[i]].contains(componentId[j]))
            {
                auto u=componentId[i],v=componentId[j];
                S[u].insert(v);
                DAG.connect(u,v);
            }
            return std::make_pair(DAG,componentId);
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

        // Decomposition of strongly connected components
        struct SCDec {
            int components; // Number of components
            std::unordered_map<UnorderedSet,int> componentId;// A map that gives the id of the strongly connected component of a node
            std::vector<UnorderedSet> topologicalOrder; // A topological order of the components
        };

        void assign(const UnorderedSet &u, std::unordered_set<UnorderedSet> &assigned,std::unordered_map<UnorderedSet,int> & cmpId,int id)
        {
            if (assigned.contains(u)) return;
            cmpId[u]=id;
            assigned.insert(u);
            for (const auto& v : reverseList[u]) assign(v,assigned,cmpId,id);
        }

        SCDec getSCComponents()
        {
            std::set<UnorderedSet> hasId;
            std::vector<UnorderedSet> topoOrder = topologicalSort();
            std::map<UnorderedSet,int> cmpId;
            int components=0;
            for (const auto& u :topoOrder) if (!hasId.contains(u))
                assign(u,hasId,cmpId,components++);
            return {components,cmpId, topoOrder};
        }

        std::pair<Graph,std::unordered_map<UnorderedSet,int>> condensationGraph()
        {
            auto [components,componentId,topologicalOrder]=getSCComponents();
            Graph DAG(components.size());
            std::vector<std::unordered_set<int>> S(components);
            for(int i=0;i<states.size();i++) for(auto j:adjacencyList[i]) if(componentId[i]!=componentId[j] && !S[componentId[i]].contains(componentId[j]))
            {
                auto u=componentId[i],v=componentId[j];
                S[u].insert(v);
                DAG.connect(u,v);
            }
            return std::make_pair(DAG,S);
        }

    };


/**
 * @brief This the class of directed Graphs
 * @details Each Graph is a couple G=(V,E) where V={0,...,n-1} are nodes, and E is a subset of VxV
 * @param n the size of the Graph
 * */
    template<typename Weight>
    struct WeightedGraph : public AbstractWeightedGraph<int,Weight>
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
                for(auto [s,_]:adjacencyList[r])
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

        int size() const override
        {
            return n;
        }

        view_or_value<AdjacentType> adjacentNodes(const int&u, bool direction) const override
        {
            if(direction)
                return std::span<const AdjacentType>(adjacencyList[u].data(),adjacencyList[u].size());
            return std::span<const AdjacentType>(reverseList[u].data(),reverseList[u].size());
        }

        view_or_value<AdjacentType> adjacentNodes(const int&u) const override
        {
            return adjacentNodes(u,true);
        }

        view_or_value<int> nodes() const override
        {
            std::vector<int> A(n);
            for(int i=0;i<n;i++)
                A[i]=i;
            return A;
        }
    };

}

#endif //CPLIBRARY_GRAPH_H
