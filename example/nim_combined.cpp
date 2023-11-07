//
// Created by ramizouari on 24/10/23.
//

#ifndef CPLIBRARY_TREE_H
#define CPLIBRARY_TREE_H

#include <optional>
#include <stdexcept>
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
//
// Created by ramizouari on 21/11/22.
//

#ifndef CPLIBRARY_UNION_FIND_H
#define CPLIBRARY_UNION_FIND_H
#include <vector>
#include <map>
#include <unordered_map>

class UnionFind
{
    int n;
    std::vector<int> parent,rank;
public:
    UnionFind(int n):n(n),rank(n),parent(n)
    {
        for(int i=0;i<n;i++)
            parent[i]=i;
    }

    void connect(int a,int b)
    {
        auto u= representative(a),v= representative(b);
        if(u==v)
            return;
        if(rank[u]<rank[u])
            parent[u]=parent[v];
        else if(rank[v]<rank[u])
            parent[v]=parent[u];
        else
        {
            parent[u]=parent[v];
            rank[v]++;
        }
    }

    int representative(int a)
    {
        if(parent[a]==a)
            return a;
        else return parent[a]= representative(parent[a]);
    }

    bool equivalent(int a,int b)
    {
        return representative(a)== representative(b);
    }
};

template<typename OrderedSet>
class OrderedUnionFind
{
    struct ParentData
    {
        const OrderedSet* parent;
        int rank;
    };
    std::map<OrderedSet,ParentData> data{};
public:
    void connect(const OrderedSet& A,const OrderedSet&B)
    {
        if(equivalent(A,B))
            return;
        auto &[C,NodeData1]=get(representative(A));
        auto &[D,NodeData2]=get(representative(B));
        auto &[p1,r1]=NodeData1;
        auto &[p2,r2]=NodeData2;
        if(r1<r2)
            p1=p2;
        else if(r1>r2)
            p2=p1;
        else
        {
            p1=p2;
            r2++;
        }
    }
    std::pair<const OrderedSet,ParentData>& get(const OrderedSet &A)
    {
        if(!data.contains(A))
        {
            auto [it,_]=data.emplace(A,ParentData{nullptr,0});
            auto &[B,NodeData]=*it;
            NodeData.parent=&B;
            return *it;
        }
        return *data.find(A);
    }
    const OrderedSet& representative(const OrderedSet&A)
    {
        auto &[B,NodeData]=get(A);
        auto &C=NodeData.parent;
        if(&B==C)
            return B;
        else {
            NodeData.parent = &representative(*NodeData.parent);
            return *NodeData.parent;
        }
    }

    bool equivalent(const OrderedSet&A,const OrderedSet&B)
    {
        return &representative(get(A).first)== &representative(get(B).first);
    }

};

template<typename UnorderedSet>
class UnorderedUnionFind
{
    struct ParentData
    {
        const UnorderedSet* parent;
        int rank;
    };
    std::unordered_map<UnorderedSet,ParentData> data{};
public:
    void connect(const UnorderedSet& A,const UnorderedSet&B)
    {
        if(equivalent(A,B))
            return;
        auto &[C,NodeData1]=get(representative(A));
        auto &[D,NodeData2]=get(representative(B));
        auto &[p1,r1]=NodeData1;
        auto &[p2,r2]=NodeData2;
        if(r1<r2)
            p1=p2;
        else if(r1>r2)
            p2=p1;
        else
        {
            p1=p2;
            r2++;
        }
    }
    std::pair<const UnorderedSet,ParentData>& get(const UnorderedSet &A)
    {
        if(!data.contains(A))
        {
            auto [it,_]=data.emplace(A,ParentData{nullptr,0});
            auto &[B,NodeData]=*it;
            NodeData.parent=&B;
            return *it;
        }
        return *data.find(A);
    }
    const UnorderedSet& representative(const UnorderedSet&A)
    {
        auto &[B,NodeData]=get(A);
        auto &C=NodeData.parent;
        if(&B==C)
            return B;
        else {
            NodeData.parent = &representative(*NodeData.parent);
            return *NodeData.parent;
        }
    }

    bool equivalent(const UnorderedSet&A,const UnorderedSet&B)
    {
        return &representative(get(A).first)== &representative(get(B).first);
    }

};



#endif //CPLIBRARY_UNION_FIND_H
//
// Created by ramizouari on 27/10/23.
//

#ifndef CPLIBRARY_GENERAL_H
#define CPLIBRARY_GENERAL_H
#include <span>
#include <variant>
#include <vector>

template<typename T>
struct view_or_value
{
    std::variant<std::vector<T>,std::span<const T>> value;
    view_or_value(std::vector<T> _value):value(std::move(_value)){}
    view_or_value(std::span<const T> _value):value(std::move(_value)){}
    std::span<const T> get() const
    {
        if(value.index()==0) return std::span<const T>(std::get<0>(value).data(),std::get<0>(value).size());
        return std::get<1>(value);
    }
    operator std::span<const T>() const
    {
        return get();
    }

    auto begin() const
    {
        return get().begin();
    }
    auto end() const
    {
        return get().end();
    }

};



namespace graph
{
    template<typename T>
    struct AbstractGraph
    {
        [[nodiscard]] virtual int size() const = 0;
        virtual view_or_value<T> adjacentNodes(const T&u, bool direction) const =0;
        virtual view_or_value<T> adjacentNodes(const T&u) const =0;
        virtual view_or_value<T> nodes() const =0;
    };

    template<typename T,typename Weight>
    struct AbstractWeightedGraph
    {
        [[nodiscard]] virtual int size() const = 0;
        using AdjacentType=std::pair<T,Weight>;
        virtual view_or_value<AdjacentType> adjacentNodes(const T&u, bool direction) const =0;
        virtual view_or_value<AdjacentType> adjacentNodes(const T&u) const =0;
        virtual view_or_value<T> nodes() const =0;
    };
}

#endif //CPLIBRARY_GENERAL_H

/**
 * @brief This the class of directed Graphs
 * @details Each Graph is a couple G=(V,E) where V={0,...,n-1} are nodes, and E is a subset of VxV
 * @param n the size of the Graph
 * */

namespace graph
{
    struct Graph : public AbstractGraph<int>
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
            for(int i=0;i<n;i++) for(auto j:adjacencyList[i]) if(!classes.equivalent(i,j) && !S[componentId[i]].contains(componentId[j]))
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
//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_DYN_SPARSE_ARRAY_H
#define CPLIBRARY_DYN_SPARSE_ARRAY_H
#include <vector>
//
// Created by ASUS on 01/12/2021.
//
#ifndef __OPERATION_H__
#define __OPERATION_H__
#include <numeric>
//
// Created by ramizouari on 01/12/2021.
//

#ifndef ACPC_PREPARATION_ABSTRACT_ALGEBRA_H
#define ACPC_PREPARATION_ABSTRACT_ALGEBRA_H
#include <complex>
#include <functional>
#include <cstdint>
using natural = std::uint64_t;
using integer = std::int64_t;
using real = long double;
using IR=real;
using IC= std::complex<IR>;
constexpr real epsilon=1e-7;

template<typename R>
R commutator(R a,R b)
{
    return a*b-b*a;
}

template<typename M,typename G=typename M::base_field>
M conj(M a)
{
    if constexpr (std::is_same_v<G, IC>)
    {
        if constexpr (std::is_same_v<G, M>)
            return std::conj(a);
        else for (auto& s : a)
            s = conj<typename std::remove_reference<decltype(s)>::type, G>(s);
    }
    return a;
}

template<typename R,typename ...StructureMetaData>
R pow(R a, long long n,StructureMetaData ... meta_info)
{
    if(n==0)
        return R(1,meta_info...);
    else if(n==1)
        return a;
    auto s=pow(a,n/2);
    return n%2?s*s*a:s*s;
}

template<typename R>
bool is_zero(const R&a)
{
    return a==R{};
}

inline bool is_zero(const IC&a)
{
    return std::abs(a) < epsilon;
}

inline bool is_zero(const real &a)
{
    return std::abs(a) < epsilon;
}

template<typename R>
R gcd(R a,R b)
{
    if(a<b)
        std::swap(a,b);
    R q,tmp;
    while(!is_zero(b))
    {
        q=a/b;
        tmp=b;
        b=a-b*q;
        a=tmp;
    }
    return a;
}

template<typename R>
R lcm(const R &a,const R &b)
{
    return a*b/gcd(a,b);
}

template<typename R=integer>
struct egcd_t
{
    R a,b,d;
};

template<typename R>
egcd_t<R> egcd(R a,R b)
{
    if(a<b)
    {
        auto e = egcd(b, a);
        std::swap(e.a,e.b);
        return e;
    }
    R q,s1=1,s2=0,t1=0,t2=1,tmp;
    while(!is_zero(b))
    {
        q=a/b;
        tmp=s2;
        s2=s1-q*s2;
        s1=tmp;
        tmp=t2;
        t2=t1-q*t2;
        t1=tmp;
        tmp=b;
        b=a-b*q;
        a=tmp;
    }
    return {s1,t1,a};
}

template<typename R>
std::pair<R,R> bezout(R a, R b)
{
    auto [u,v,_]=egcd(a,b);
    return {u,v};
}

template<typename B>
B next_gray(B n)
{
    return n^(n>>1);
}

template<typename F,typename R>
std::pair<integer,integer> floyd_functional_cycle(F && f,R x0)
{
    /*
     * Find a period
     * */
    R x=x0,y=x;
    integer m=0;
    do
    {
        x=f(x);
        y=f(f(y));
        m++;
    }while(y!=x);
    /*
     * Find offset
     * */
    x=x0,y=x;
    for(int i=0;i<m;i++)
        y=f(y);
    int offset=0;
    while(x!=y)
    {
        x=f(x);
        y=f(y);
        offset++;
    }

    /*
     * Find fundamental period
     * */
    y=f(x);
    integer period=1;
    while(x!=y) {
        y = f(y);
        period++;
    }
    return std::make_pair(period,offset);
}


template<typename F,typename R>
integer functional_period(F &&f, R x)
{
    /*
    * Find a period
    * */
    R y=x;
    integer m=0;
    do
    {
        x=f(x);
        y=f(f(y));
        m++;
    }while(y!=x);
    return m;
}


#endif //ACPC_PREPARATION_ABSTRACT_ALGEBRA_H
#include <memory>

template<typename T>
struct binary_operation
{
    using type=T;
    template<typename H0,typename ...H>
    T operator()(const H0&a,const H&... b) const
    {
        if constexpr (sizeof...(b) == 0)
            return a;
        else return reduce(a,this->operator()(b...));
    }
    virtual T reduce(const T& a, const T& b) const = 0;
    virtual T neutral_element() const
    {
        return T{};
    }
};

template<typename T>
struct invertible_operation
{
    virtual T inv(const T& a) const = 0;
};

template<typename T>
struct monoid_plus_t:public binary_operation<T> {
    T reduce(const T &a, const T &b) const override {
        return a + b;
    }
    inline static T neutral{};
};

template<typename T>
struct plus_t:public monoid_plus_t<T>,public invertible_operation<T>
{
    T inv(const T&a) const override
    {
        return -a;
    }
};

template<typename T>
struct multiplies_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return a*b;
    }

    inline static T neutral=T(1);
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct field_multiplies_t:public multiplies_t<T>,public invertible_operation<T>
{
    T inv(const T&a) const
    {
        return a.inv();
    }

};

template<>
struct field_multiplies_t<real>:public multiplies_t<real>,public invertible_operation<real>
{
    real inv(const real& a)const
    {
        return 1./a;
    }
};

template<>
struct field_multiplies_t<IC>:public multiplies_t<IC>,public invertible_operation<IC>
{
    IC inv(const IC& a)const
    {
        return IC(1)/a;
    }
};

template<typename T>
struct max_t:public binary_operation<T>
{
    T e;
    explicit max_t(T _e):e(_e){}
    max_t(): max_t(T{}){}
    T reduce(const T&a,const T&b) const override
    {
        return std::max(a,b);
    }

    inline static T neutral{0};
    T neutral_element() const override
    {
        return e;
    }
};

template<typename T>
struct min_t:public binary_operation<T>
{
    T e;
    explicit min_t(T _e):e(_e){}
    min_t(): min_t(T{}){}

    T reduce(const T&a,const T&b) const override
    {
        return std::min(a,b);
    }

    inline static T neutral{};

    T neutral_element() const override
    {
        return e;
    }
};

template<typename T>
struct gcd_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return gcd(a,b);
    }

    inline static T neutral{0};
};

template<typename T>
struct lcm_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return lcm(a,b);
    }

    inline static T neutral{1};
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct xor_t:public binary_operation<T>,public invertible_operation<T>
{
    T reduce(const T&a,const T&b) const
    {
        return a^b;
    }

    T inv(const T&a) const
    {
        return a;
    }

    inline static T neutral{};
};

template<typename T>
struct and_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return a&b;
    }

    inline static T neutral=static_cast<T>(-1);
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct or_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return a|b;
    }

    inline static T neutral{};
};

template<typename T>
struct logical_and_t :public binary_operation<T>
{
    T reduce(const T& a, const T& b) const override
    {
        return a && b;
    }

    inline static T neutral{true};
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct logical_or_t :public binary_operation<T>
{
    T reduce(const T& a, const T& b) const override
    {
        return a || b;
    }

    inline static T neutral{false};
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct logical_xor_t :public binary_operation<T>,public invertible_operation<T>
{
    T reduce(const T& a, const T& b) const override
    {
        return !a && b || a && !b;
    }
    T inv(const T&a) const
    {
        return !a;
    }
    inline static T neutral{false};
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
class binary_operation_ptr
{
    std::shared_ptr<binary_operation<T>> op;
public:
    binary_operation_ptr(std::shared_ptr<binary_operation<T>> value): op(value){}
    template<typename ...H>
    auto operator()(const H&... h) const
    {
        return op->operator()(h...);
    }

    auto neutral_element() const
    {
        return op->neutral_element();
    }
};

template<typename T>
class invertible_binary_operation_ptr : public binary_operation_ptr<T>
{
    std::shared_ptr<invertible_operation<T>> inverter;
public:
    invertible_binary_operation_ptr(std::shared_ptr<binary_operation<T>> b,
                                    std::shared_ptr<invertible_operation<T>> I): binary_operation_ptr<T>(b),inverter(I){}
    using binary_operation_ptr<T>::operator();
    using binary_operation_ptr<T>::neutral_element;
    auto inv(const T& a) const
    {
        return inverter->inv(a);
    }
};

#endif
#include <memory>
//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_BITS_H
#define CPLIBRARY_BITS_H

inline unsigned int bit_log(unsigned int n)
{
    unsigned char a=0,b=30,r=0;
    while(a<=b)
    {
        auto c=(a+b)/2;
        if(n>>c)
            a=c+1;
        else
        {
            b=c-1;
            r=c-1;
        }
    }
    if(r && (1<<(r-1))==n)
        return r-1;
    return r;
}

inline unsigned int bit_floor(unsigned int n)
{
    return 1<<bit_log(n);
}

inline unsigned int bit_ceil(unsigned int n)
{
    return 1<<(bit_log(n)+1);
}

#endif //CPLIBRARY_BITS_H

namespace data_structures::dynamic
{
    template<typename T>
    struct sparse_array
    {
        int n,h;
        std::vector<std::vector<T>> S;
        binary_operation_ptr<T> F;
    public:
        sparse_array(const std::vector<T>&A, std::shared_ptr<binary_operation<T>> _F):n(bit_ceil(A.size())),h(bit_log(n)),S(h+1),F(_F)
        {
            int r=1;
            for(int i=h;i>=0;i--,r*=2)
                S[i].resize(n-r+1,F->neutral_element());
            for(int i=0;i<A.size();i++)
                S[h][i]=A[i];
            r=1;
            for(int i=h-1;i>=0;i--,r*=2) for(int j=0;j<=n-2*r;j++)
                    S[i][j]=F(S[i+1][j],S[i+1][j+r]);
        }

        T query(int l,int r) const
        {
            if(l>=r)
                return F->neutral_element();
            auto d=r-l;
            auto s=bit_floor(d);
            auto b=bit_log(s);
            return F(S[h-b][l],S[h-b][r-s]);
        }
    };
}

#endif //CPLIBRARY_SPARSE_ARRAY_H
//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_FIXED_SPARSE_ARRAY_H
#define CPLIBRARY_FIXED_SPARSE_ARRAY_H
#include <vector>
namespace data_structures::fixed
{
    template<typename O>
    struct sparse_array
    {
        using T=typename O::type;
        using type = T;
        inline static O F=O();
        int n,h;
        std::vector<std::vector<T>> S;
    public:
        sparse_array(const std::vector<T>&A):n(bit_ceil(A.size())),h(bit_log(n)),S(h+1)
        {
            int r=1;
            for(int i=h;i>=0;i--,r*=2)
                S[i].resize(n-r+1,O::neutral);
            for(int i=0;i<A.size();i++)
                S[h][i]=A[i];
            r=1;
            for(int i=h-1;i>=0;i--,r*=2) for(int j=0;j<=n-2*r;j++)
                    S[i][j]=F(S[i+1][j],S[i+1][j+r]);
        }

        T query(int l,int r) const
        {
            if(l>=r)
                return O::neutral;
            auto d=r-l;
            auto s=bit_floor(d);
            auto b=bit_log(s);
            return F(S[h-b][l],S[h-b][r-s]);
        }
    };


}
#endif //CPLIBRARY_SPARSE_ARRAY_H
#include <memory>
#include <queue>
#include <random>
#include <chrono>

namespace graph
{

    namespace ds= data_structures;

    struct HLDIndex
    {
        int hld_id;
        int index;
        HLDIndex(int _hld_id, int _index): hld_id(_hld_id), index(_index){}
    };

    template<typename Weight>
    struct HeavyLightDecomposition
    {
        std::vector<bool> is_heavy;
        std::vector<std::pair<int,int>> heavy_path_endpoints;
        std::vector<int> component_size;
        std::vector<HLDIndex> HLD_mapper;
        std::vector<std::vector<Weight>> components;
    };

    template<>
    struct HeavyLightDecomposition<void>
    {
        std::vector<bool> is_heavy;
        std::vector<std::pair<int,int>> heavy_path_endpoints;
        std::vector<int> component_size;
        std::vector<HLDIndex> HLD_mapper;
    };

    enum class TreeStats
    {
        NONE=0,
        SIZE=0b00001,
        HEAVY_EDGES=0b00011,
        LCA=0b00100,
        HLD=0b01111,
        RANGE_QUERIES=0b11111,
    };
    bool operator &(TreeStats a,TreeStats b)
    {
        return static_cast<int>(a) & static_cast<int>(b);
    }

    struct Tree : Graph
    {

        bool reversed=false;
        std::vector<int> subtree_size;
        std::vector<std::optional<int>> parent;
        int root;
        HeavyLightDecomposition<void> HLD;
        Tree(int n,int _root):Graph(n),root(_root),subtree_size(n),parent(n)
        {
            HLD.is_heavy.resize(n,false);
        }
        explicit Tree(int n):Tree(n,0){}

        void setParent(int u,int v)
        {
            if(reversed) std::swap(u,v);
            this->connect(u,v);
            parent[u].emplace(v);
        }

        std::vector<int> &children(int u)
        {
            return reversed?this->adjacencyList[u] : this->reverseList[u] ;
        }

        [[nodiscard]] const std::vector<int> &children(int u) const
        {
            return reversed?this->adjacencyList[u] : this->reverseList[u] ;
        }

        void buildStatistics(TreeStats stats=TreeStats::HLD)
        {
            updateRoot();
            if(stats & TreeStats::SIZE) updateSize(root);
            //TODO: Optimize heavy edges
            if(stats & TreeStats::HEAVY_EDGES) updateHeavyEdges(root);
            if(stats & TreeStats::LCA) buildLCA();
            auto t1=std::chrono::high_resolution_clock::now();
            //TODO: Optimize HLD
            if(stats & TreeStats::HLD) buildHeavyLightDecomposition();
        }

        void adjacentReRoot(int new_root)
        {
            if(parent[new_root] != root)
                throw std::invalid_argument("new root must be adjacent to old root");
            auto u=*parent[new_root];
            parent[new_root]=std::nullopt;
            parent[root].emplace(new_root);
            auto delta=subtree_size[new_root];
            subtree_size[new_root]=subtree_size[root];
            subtree_size[root]-=delta;
            root=new_root;
        }

        void reRoot(int new_root)
        {
            std::queue<int> Q;
            std::vector<bool> visited(n);
            Q.push(new_root);
            visited[new_root]=true;
            std::vector<std::vector<int>> newAdjacencyList(n),newReverseList(n);
            while(!Q.empty())
            {
                auto u=Q.front();
                Q.pop();
                for(auto v:this->adjacencyList[u]) if(!visited[v])
                    {
                        visited[v]=true;
                        Q.push(v);
                        newReverseList[u].emplace_back(v);
                        newAdjacencyList[v].emplace_back(u);
                    }
                for(auto v:this->reverseList[u]) if(!visited[v])
                    {
                        visited[v]=true;
                        Q.push(v);
                        newReverseList[u].emplace_back(v);
                        newAdjacencyList[v].emplace_back(u);
                    }
            }
            this->adjacencyList=std::move(newAdjacencyList);
            this->reverseList=std::move(newReverseList);
            updateRoot();
        }

        int leastCommonAncestor(int u,int v)
        {
            if(lca_data)
            {
                auto [a,b]=euler_tour_endpoints[u];
                auto [c,d]=euler_tour_endpoints[v];
                if(a>c)
                {
                    std::swap(a,c);
                    std::swap(b,d);
                }
                if(b<d)
                    return lca_data->query(a,c).second;
                else
                    return lca_data->query(a,d).second;
            }
            else
            {
                while(u!=v)
                {
                    if(subtree_size[u]>subtree_size[v])
                        u=*parent[u];
                    else
                        v=*parent[v];
                }
                return u;
            }
        }
        void updateSize(int u)
        {
            subtree_size[u]=1;
            for(auto v:children(u))
            {
                updateSize(v);
                subtree_size[u]+=subtree_size[v];
            }
        }

        void updateRoot()
        {
            for(int i=0;i<n;i++)
                parent[i]=this->adjacencyList[i].empty()?std::nullopt:std::make_optional(this->adjacencyList[i][0]);
            for(int i=0;i< n;i++)
                if(!parent[i])
                    root=i;
        }

        void updateHeavyEdges(int u)
        {
            for(int i=0;i<children(u).size();i++)
            {
                auto v=children(u)[i];
                if(subtree_size[v]>=(subtree_size[u]+1)/2)
                    HLD.is_heavy[v]=true;
                updateHeavyEdges(v);
            }
        }
        void buildLCA()
        {
            std::vector<HeightData> A;
            euler_tour_endpoints.resize(n);
            eulerTour(root,0,A);
            min_t<HeightData>::neutral.first=std::numeric_limits<int>::max();
            lca_data=std::make_unique<ds::fixed::sparse_array<min_t<HeightData>>>(A);
        }

        void buildHeavyLightDecomposition()
        {
            std::vector<int> stack;
            stack.push_back(root);
            HLD.HLD_mapper.resize(n, HLDIndex(-1, -1));
            HLD.HLD_mapper[root]={0, 0};
            HLD.heavy_path_endpoints.emplace_back(root,root);
            int components=1;
            while(!stack.empty())
            {
                auto u=stack.back();
                stack.pop_back();
                for(auto v:children(u))
                {
                    if(HLD.is_heavy[v])
                    {
                        auto &[_,y] = HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id];
                        stack.push_back(v);
                        HLD.HLD_mapper[v]={HLD.HLD_mapper[u].hld_id, HLD.HLD_mapper[u].index + 1};
                        y=v;
                    }
                    else
                    {
                        HLD.heavy_path_endpoints.emplace_back(v,v);
                        stack.push_back(v);
                        HLD.HLD_mapper[v]=HLDIndex(components, 0);
                        components++;
                    }
                }
            }
        }

        int distance(int a,int b)
        {
            auto lca=leastCommonAncestor(a,b);
            return distance_with_lca(a,lca)+distance_with_lca(b,lca);
        }

        int distance_with_lca(int u,int lca)
        {
            int d=0;
            while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id)
            {
                if(HLD.is_heavy[u])
                {
                    d+=HLD.HLD_mapper[u].index;
                    u=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id].first;
                }
                else
                {
                    d++;
                    u=*parent[u];
                }
            }
            return d+= HLD.HLD_mapper[u].index - HLD.HLD_mapper[lca].index;
        }

        int centroid()
        {
            reRoot(centroid(root,root));
            return root;
        }

    protected:
        int centroid(int u,std::optional<int> p)
        {
            for(auto v:children(u)) if(v!= p && subtree_size[v]>=(subtree_size[u]+1)/2)
                {
                    adjacentReRoot(v);
                    return centroid(v,u);
                }
            return u;
        }
        using HeightData=std::pair<int,int>;
        using EnpointsData = std::pair<int,int>;
        std::unique_ptr<ds::fixed::sparse_array<min_t<HeightData>>> lca_data;
        std::vector<EnpointsData> euler_tour_endpoints;
        void eulerTour(int u,int height,std::vector<HeightData> &A)
        {
            euler_tour_endpoints[u].first=A.size();
            for(auto v: children(u))
            {
                A.emplace_back(height,u);
                eulerTour(v,height+1,A);
            }
            A.emplace_back(height,u);
            euler_tour_endpoints[u].second=A.size();
        }
    };

    template<typename Weight>
    struct WeightedTree : public WeightedGraph<Weight>
    {

        bool reversed=false;
        std::vector<int> subtree_size;
        using AdjacentType=WeightedGraph<Weight>::AdjacentType;
        std::vector<std::optional<AdjacentType>> parent;
        int root;
        HeavyLightDecomposition<Weight> HLD;
        WeightedTree(int n,int _root):WeightedGraph<Weight>(n),root(_root),subtree_size(n),parent(n)
        {
            HLD.is_heavy.resize(n);
        }
        explicit WeightedTree(int n):WeightedTree(n,0){}

        void setParent(int u,int v,Weight w)
        {
            if(reversed) std::swap(u,v);
            this->connect(u,v,w);
            parent[u].emplace(v,w);
        }

        std::vector<AdjacentType> &children(int u)
        {
            return reversed?this->adjacencyList[u] : this->reverseList[u] ;
        }

        void buildStatistics(TreeStats stats=TreeStats::HLD)
        {
            updateRoot();
            if(stats & TreeStats::SIZE) updateSize(root);
            if(stats & TreeStats::HEAVY_EDGES) updateHeavyEdges(root);
            if(stats & TreeStats::LCA) buildLCA();
            if(stats & TreeStats::HLD) buildHeavyLightDecomposition();
        }

        void adjacentReRoot(int new_root)
        {
            if(!parent[new_root] || parent[new_root]->first != root)
                throw std::invalid_argument("new root must be adjacent to old root");
            auto [u,w]=*parent[new_root];
            parent[new_root]=std::nullopt;
            parent[root].emplace(new_root,w);
            auto delta=subtree_size[new_root];
            subtree_size[new_root]=subtree_size[root];
            subtree_size[root]-=delta;
        }

        void reRoot(int new_root)
        {
            std::queue<int> Q;
            std::vector<bool> visited(WeightedGraph<Weight>::n);
            Q.emplace(new_root);
            visited[new_root]=true;
            std::vector<std::vector<AdjacentType>> newAdjacencyList(WeightedGraph<Weight>::n),newReverseList(WeightedGraph<Weight>::n);
            while(!Q.empty())
            {
                auto u=Q.front();
                Q.pop();
                for(auto [v,w]:this->adjacencyList[u]) if(!visited[v])
                    {
                        visited[v]=true;
                        Q.emplace(v);
                        newReverseList[u].emplace_back(v,w);
                        newAdjacencyList[v].emplace_back(u,w);
                    }
                for(auto [v,w]:this->reverseList[u]) if(!visited[v])
                    {
                        visited[v]=true;
                        Q.emplace(v);
                        newReverseList[u].emplace_back(v,w);
                        newAdjacencyList[v].emplace_back(u,w);
                    }
            }
            this->adjacencyList=std::move(newAdjacencyList);
            this->reverseList=std::move(newReverseList);
            updateRoot();
        }

        int leastCommonAncestor(int u,int v)
        {
            if(lca_data)
            {
                auto [a,b]=euler_tour_endpoints[u];
                auto [c,d]=euler_tour_endpoints[v];
                if(a>c)
                {
                    std::swap(a,c);
                    std::swap(b,d);
                }
                if(b<d)
                    return lca_data->query(a,c).second;
                else
                    return lca_data->query(a,d).second;
            }
            else
            {
                while(u!=v)
                {
                    if(subtree_size[u]>subtree_size[v])
                        u=parent[u]->first;
                    else
                        v=parent[v]->first;
                }
                return u;
            }
        }
        void updateSize(int u)
        {
            subtree_size[u]=1;
            for(auto [v,_]:children(u))
            {
                updateSize(v);
                subtree_size[u]+=subtree_size[v];
            }
        }

        void updateRoot()
        {
            for(int i=0;i<WeightedGraph<Weight>::n;i++)
                parent[i]=this->adjacencyList[i].empty()?std::nullopt:std::make_optional(this->adjacencyList[i][0]);
            for(int i=0;i<WeightedGraph<Weight>::n;i++)
                if(!parent[i])
                    root=i;
        }

        void updateHeavyEdges(int u)
        {
            for(int i=0;i<children(u).size();i++)
            {
                auto [v,_]=children(u)[i];
                if(subtree_size[v]>=(subtree_size[u]+1)/2)
                    HLD.is_heavy[v]=true;
                updateHeavyEdges(v);
            }
        }
        void buildLCA()
        {
            std::vector<HeightData> A;
            euler_tour_endpoints.resize(WeightedGraph<Weight>::n);
            eulerTour(root,0,A);
            min_t<HeightData>::neutral.first=std::numeric_limits<int>::max();
            lca_data=std::make_unique<ds::fixed::sparse_array<min_t<HeightData>>>(A);
        }

        void buildHeavyLightDecomposition()
        {
            std::vector<int> stack;
            stack.push_back(root);
            HLD.HLD_mapper.resize(WeightedTree<Weight>::n, HLDIndex(-1, -1));
            HLD.HLD_mapper[root]={0, 0};
            HLD.components.emplace_back();
            HLD.heavy_path_endpoints.emplace_back(root,root);
            while(!stack.empty())
            {
                auto u=stack.back();
                stack.pop_back();
                for(auto [v,w]:children(u))
                {
                    if(HLD.is_heavy[v])
                    {
                        auto &C=HLD.components[HLD.HLD_mapper[u].hld_id];
                        auto &[_,y] = HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id];
                        stack.push_back(v);
                        HLD.HLD_mapper[v]={HLD.HLD_mapper[u].hld_id, HLD.HLD_mapper[u].index + 1};
                        y=v;
                        C.push_back(w);
                    }
                    else
                    {
                        HLD.heavy_path_endpoints.emplace_back(v,v);
                        stack.push_back(v);
                        HLD.HLD_mapper[v]=HLDIndex(HLD.components.size(), 0);
                        HLD.components.emplace_back();
                    }
                }
            }
            for(const auto &C:HLD.components)
                HLD.component_size.push_back(C.size());
        }

        int distance(int a,int b)
        {
            auto lca=leastCommonAncestor(a,b);
            return distance_with_lca(a,lca)+distance_with_lca(b,lca);
        }

        int distance_with_lca(int u,int lca)
        {
            int d=0;
            while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id)
            {
                if(WeightedTree<Weight>::HLD.is_heavy[u])
                {
                    d+=HLD.HLD_mapper[u].index;
                    u=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id].first;
                }
                else
                {
                    d++;
                    u=parent[u]->first;
                }
            }
            return d+= HLD.HLD_mapper[u].index - HLD.HLD_mapper[lca].index;
        }


    protected:
        using HeightData=std::pair<int,int>;
        using EnpointsData = std::pair<int,int>;
        std::unique_ptr<ds::fixed::sparse_array<min_t<HeightData>>> lca_data;
        std::vector<EnpointsData> euler_tour_endpoints;
        void eulerTour(int u,int height,std::vector<HeightData> &A)
        {
            euler_tour_endpoints[u].first=A.size();
            for(auto [v,_]: children(u))
            {
                A.emplace_back(height,u);
                eulerTour(v,height+1,A);
            }
            A.emplace_back(height,u);
            euler_tour_endpoints[u].second=A.size();
        }
    };

}


#endif //CPLIBRARY_TREE_H
//
// Created by ramizouari on 07/11/23.
//

#ifndef CPLIBRARY_ISOMORPHISM_H
#define CPLIBRARY_ISOMORPHISM_H

namespace graph
{
    using char_deque=std::deque<char>;

    std::shared_ptr<char_deque> string_encode(const graph::Tree & T, int u)
    {
        std::vector<std::shared_ptr<char_deque>> X;
        const auto &C=T.children(u);
        if(C.empty())
            return std::make_shared<char_deque>(char_deque{'(', ')'});
        X.reserve(C.size());
        for(auto v:C) X.push_back(string_encode(T,v));
        std::sort(X.begin(),X.end(),[](const auto &x,const auto &y)
        {
            return *x < *y;
        });
        auto it=std::max_element(X.begin(),X.end(),[](const auto &x,const auto &y) {
            return x->size() < y->size();
        });
        auto r=std::distance(X.begin(),it);
        std::shared_ptr<char_deque> Z=X[r];
        for(int i=r-1;i>=0;i--)
            std::copy(X[i]->rbegin(),X[i]->rend(),std::front_inserter(*Z));
        for(int i=r+1;i<X.size();i++)
            std::copy(X[i]->begin(),X[i]->end(),std::back_inserter(*Z));
        Z->push_back(')');
        Z->push_front('(');
        return Z;
    }

    std::string string_encode(const graph::Tree &T)
    {
        auto E=string_encode(T,T.root);
        return std::string(E->begin(),E->end());
    }

    std::optional<int> second_centroid(graph::Tree & T)
    {
        auto u=T.root;
        auto X=T.children(u);
        for(auto v:X)
        {
            bool is_centroid=true;
            T.adjacentReRoot(v);
            if(T.subtree_size[u]>T.subtree_size[v]/2)
                is_centroid=false;
            T.adjacentReRoot(u);
            if(is_centroid)
                return v;
        }
        return std::nullopt;
    }

    std::vector<std::string> full_string_encoding(graph::Tree &T)
    {
        T.centroid();
        T.buildStatistics(graph::TreeStats::SIZE);
        std::vector<std::string> A;
        A.push_back(string_encode(T));
        auto p= second_centroid(T);
        if(p.has_value())
        {
            T.reRoot(*p);
            T.buildStatistics(graph::TreeStats::SIZE);
            A.push_back(string_encode(T));
            if(A[0] > A[1])
                std::swap(A[0],A[1]);
        }
        return A;
    }

    struct TreeHolder
    {
        std::shared_ptr<Tree> tree;
        TreeHolder(Tree&& t): tree(std::make_shared<Tree>(std::move(t))){}
        TreeHolder(int n): tree(std::make_shared<Tree>(n)){}

        Tree* operator->() const
        {
            return tree.get();
        }
        Tree& operator*() const
        {
            return *tree;
        }
    };


    struct IsoTreeEq
    {
        bool operator()(Tree &a, Tree &b) const
        {
            if(a.size() != b.size())
                return false;
            a.centroid();
            b.centroid();
            auto M1= string_encode(a);
            auto M2= string_encode(b);
            if(M1==M2)
                return true;
            auto c= second_centroid(a);
            if(c.has_value())
            {
                a.reRoot(*c);
                auto M3= string_encode(a);
                return M3==M2;
            }
            return false;
        }

        bool operator()(const TreeHolder &a, const TreeHolder &b) const
        {
            return (*this)(*a,*b);
        }
    };

    struct IsoTreeCmp
    {
        bool operator()(Tree &a, Tree &b) const
        {
            if(a.size()!=b.size())
                return a.size() < b.size();
            auto X=full_string_encoding(a);
            auto Y=full_string_encoding(b);
            for(int i:{0,1}) for(int j:{0,1}) if(i<X.size() && j<Y.size() && X[i]==Y[j])
            {
                if(i!=0 || j!=0)
                    throw std::runtime_error("?????????????????????????????????????");
                else return false;
            }
            return X < Y;
        }

        bool operator()(const TreeHolder& a, const TreeHolder& b) const
        {
            return (*this)(*a,*b);
        }
    };

    struct IsoTreeHash
    {
        std::hash<std::string> H;
    public:
        size_t operator()(Tree &a) const
        {
            auto A= full_string_encoding(a);
            return std::accumulate(A.begin(),A.end(),0ULL,[&H=this->H](auto x,auto y){
                return x^H(y);
            });
        }

        size_t operator()(const TreeHolder& a) const
        {
            return (*this)(*a);
        }
    };
}

#endif //CPLIBRARY_ISOMORPHISM_H
#include <iostream>

int main()
{
    std::ios_base::sync_with_stdio(false);
    int T;
    std::cin >> T;
    graph::IsoTreeEq treeEq;
    using TreesSet=std::set<graph::TreeHolder,graph::IsoTreeCmp>;
    for(int t=1;t<=T;t++)
    {
        TreesSet Z;
        int n;
        std::cin >> n;
        graph::TreeHolder T1(n), T2(n);
        for(int i=0;i<n-1;i++)
        {
            int u,v;
            std::cin >> u >> v;
            u--;
            v--;
            T1->connect(u,v);
        }
        for(int i=0;i<n-1;i++)
        {
            int u,v;
            std::cin >> u >> v;
            u--;
            v--;
            T2->connect(u,v);
        }
        T1->reRoot(0);
        T2->reRoot(0);
        T1->buildStatistics(graph::TreeStats::SIZE);
        T2->buildStatistics(graph::TreeStats::SIZE);
        Z.emplace(std::move(T1));
        Z.emplace(std::move(T2));
        auto result=(Z.size()==1?"YES":"NO");
        std::cout << result << '\n';
    }
}
