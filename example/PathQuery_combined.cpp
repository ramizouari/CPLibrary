#include <iostream>
//
// Created by ramizouari on 25/10/23.
//

#ifndef CPLIBRARY_RANGE_QUERIES_H
#define CPLIBRARY_RANGE_QUERIES_H
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
//
// Created by ASUS on 01/12/2021.
//
#ifndef __DATA_STRUCTURES_H__
#define __DATA_STRUCTURES_H__
#include <vector>
//
// Created by ramizouari on 17/11/2021.
//

#ifndef __STATISTIC_NODE__
#define __STATISTIC_NODE__
#include <tuple>
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

#endif
#include <variant>

/**
* @definition Order Statistic Tree: It is an AVL-tree augmented with a statistic.
* The statistic is a function of:
 * </ol>
* <li> The key
* <li> The value
* <li> left and right subtree
 * </ol>
* @Requirements
* <ol><li> a class T
* <li> (T,<=) is a totally ordered set
* <li> S is a statistic type
* </ol>
* @tparam T the type of the key
* @tparam V the type of the value
* @tparam S the statistic type: a statistic S is a class having some attributes serving as additional informations
* Those informations can be aggregated via the static update method
* @Requirements
* <ol>
* <li> class S
* <li> S has a public static method called update which accepts Tree<T,V,S>.
* <li> the update method "updates" adequately the statistics, and only the statistics
 * </ol>
*/
template<typename T,typename V,typename S>
struct statistic_node;

template<typename T,typename V,typename S>
struct statistic_node
{
    T v;
    V data;
    int h;
    S statistic;
    statistic_node *left,*right,*parent;
    statistic_node(T _v,V _data,statistic_node* _parent=nullptr):v(_v),data(_data),left(nullptr),right(nullptr),parent(_parent),h(1),statistic(v,data){}
    void update()
    {
        h=std::max(left?left->h:0,right?right->h:0)+1;
        S::update(this);
    }
};


/*
* Get the height of the (sub)-tree
*/
template<typename T,typename V,typename S>
int height(statistic_node<T,V,S>*node)
{
    return node?node->h:0;
}


/*
* Get the balance of the current node
*/
template<typename T,typename V,typename S>
int balance(statistic_node<T, V, S>* tree)
{
    return height(tree->left) - height(tree->right);
}

/**
* Find node that has the given key, return nullptr otherwise
* @Notes
* 1. If the value is altered, and the statistic depends on that value. The calculated statistic will not be updated.
*   In that case, better use the insert_or_assign function or the update function
* 2. If the key is altered, It will probably cause the tree to be in an inconsistent state unless the tree does not have a key
* on the interval whose limit points are the old key value, and the new one.
*/
template<typename T, typename V, typename S>
statistic_node<T,V,S>* find(statistic_node<T, V, S>* node, const typename std::common_type<T>::type& v)
{
    if (!node)
        return nullptr;
    if (node->v == v)
        return node->data;
    else if (node->v < v)
        return find(node->right, v);
    else return find(node->left, v);
}

/*
* Find the data mapped by the given key
* @Requirements
* The Order Statistic Tree must contain at least one such key
* @Exception
* std::logic_error thrown if no such key is found
* @Notes
* If the value is altered, and the statistic depends on that value. The calculated statistic will not be updated.
* In that case, better use the insert_or_assign function or the update function
*/
template<typename T, typename V, typename S>
V& value_at(statistic_node<T, V, S>* node, const typename std::common_type<T>::type& v)
{
    [[unlikely]]
    if (!node)
        throw std::out_of_range("key does not exist");
    if (node->v == v)
        return node;
    else if (node->v < v)
        return value_at(node->right, v);
    else return value_at(node->left, v);
}

/*
* Get the node whose key is strictly bigger than v
*/
template<typename T,typename V,typename S>
statistic_node<T,V,S>* upper_bound(statistic_node<T,V,S>* tree,const 
    typename std::common_type<T>::type& v)
{
    if(tree==nullptr)
        return nullptr;
    if(tree->v<=v)
        return upper_bound(tree->right,v);
    else
    {
        auto w=upper_bound(tree->left,v);
        if(w==nullptr)
            return tree;
        return w;
    }
}

/*
* Get the node whose key is strictly smaller than v
*/
template<typename T,typename V,typename S>
statistic_node<T,V,S>* reverse_upper_bound(statistic_node<T,V,S>* tree,
    const typename std::common_type<T>::type & v)
{
    if(tree==nullptr)
        return nullptr;
    if(tree->v>=v)
        return reverse_upper_bound(tree->left,v);
    else
    {
        auto w=reverse_upper_bound(tree->right,v);
        if(w==nullptr)
            return tree;
        return w;
    }
}

/*
* Get the node whose key is not smaller than v
*/
template<typename T,typename V,typename S>
statistic_node<T,V,S>* lower_bound(statistic_node<T,V,S>* tree,
    const typename std::common_type<T>::type& v)
{
    if(tree==nullptr)
        return nullptr;
    if(tree->v<v)
        return lower_bound(tree->right,v);
    else
    {
        auto w=lower_bound(tree->left,v);
        if(w==nullptr)
            return tree;
        return w;
    }
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* reverse_lower_bound(statistic_node<T,V,S>* tree,
    const typename std::common_type<T>::type &v)
{
    if(tree==nullptr)
        return nullptr;
    if(tree->v>v)
        return reverse_lower_bound(tree->left,v);
    else
    {
        auto w=reverse_lower_bound(tree->right,v);
        if(w==nullptr)
            return tree;
        return w;
    }
}

/*
* Get the node with the smallest key
*/
template<typename T, typename V, typename S>
statistic_node<T, V, S>* begin(statistic_node<T, V, S>* tree)
{
    if (!tree)
        return nullptr;
    while (tree->left)
        tree = tree->left;
    return tree;
}

/*
* Get the successor of the current node
*/
template<typename T, typename V, typename S>
statistic_node<T, V, S>* next(statistic_node<T, V, S>* tree)
{
    if (tree == nullptr)
        return nullptr;
    if (tree->right)
    {
        tree = tree->right;
        while (tree->left)
            tree = tree->left;
        return tree;
    }
    else
    {
        auto tmp = tree;
        tree = tree->parent;
        while (tree && tree->v < tmp->v)
            tree = tree->parent;
        if (!tree)
            return nullptr;
        return tree;
    }
}

/*
* Get the previous of the current node
*/
template<typename T, typename V, typename S>
statistic_node<T, V, S>* prev(statistic_node<T, V, S>* tree)
{
    if (tree == nullptr)
        return nullptr;
    if (tree->left)
    {
        tree = tree->left;
        while (tree->right)
            tree = tree->right;
        return tree;
    }
    else
    {
        auto tmp = tree;
        tree = tree->parent;
        while (tree && tree->v > tmp->v)
            tree = tree->parent;
        if (!tree)
            return nullptr;
        return tree;
    }
}

/*
* Applies a right rotation to the ordered statistic tree on the current node.
*/
template<typename T,typename V,typename S>
statistic_node<T,V,S>* rebalance_right(statistic_node<T,V,S>* x)
{
    auto y=x->left,B=y->right;
    y->right=x;
    y->parent=x->parent;
    if(x->parent)
    {
        if(x->parent->left==x)
            x->parent->left=y;
        else x->parent->right=y;
    }
    x->parent=y;
    x->left=B;
    if(B) B->parent=x;
    x->update();
    return y;
}

/*
* Applies a left rotation to the ordered statistic tree on the current node.
*/
template<typename T,typename V,typename S>
statistic_node<T,V,S>* rebalance_left(statistic_node<T,V,S>* x)
{
    auto y=x->right,B=y->left;
    y->left=x;
    y->parent=x->parent;
    if(x->parent)
    {
        if (x->parent->left == x)
            x->parent->left = y;
        else x->parent->right = y;
    }
    x->parent=y;
    x->right=B;
    if(B) B->parent=x;
    x->update();
    return y;
}

/*
* Rebalance the ordered statistic tree. 
*/
template<typename T,typename V,typename S>
statistic_node<T,V,S>* rebalance(statistic_node<T,V,S>* x)
{
    if(!x)
        return nullptr;
    if(balance(x)<-1)
    {
        if(balance(x->right)==1)
            rebalance_right(x->right);
        x=rebalance_left(x);
    }
    else if(balance(x)>1)
    {
        if(balance(x->left)==-1)
            rebalance_left(x->left);
        x=rebalance_right(x);
    }
    x->update();
    if(!x->parent)
        return x;
    return rebalance(x->parent);
}

/*
* Insert (v,data) into the ordered statistic tree
* @Cases
* 1. If or_assign is false, (v,data) is inserted while allowing duplicates
* 2. if or_assign is true, (v,data) is inserted if the tree does not have a key v. Otherwise the value mapped by v is changed to data.
*/
template<typename T,typename V,typename S>
statistic_node<T,V,S>* insert(statistic_node<T,V,S>* tree,const typename std::common_type<T>::type& v,
    const typename std::common_type<V>::type& data,bool or_assign=false)
{
    if(!tree)
    {
        tree = new statistic_node<T,V,S>(v,data);
        return tree;
    }
    auto p=lower_bound(tree,v);
    if(p==nullptr)
    {
        p=tree;
        while(p->right) p=p->right;
    }
    else if (or_assign && p->v == v)
    {
        p->data = data;
        p->update();
        return rebalance(p);
    }
    else if(p->left)
    {
        p=p->left;
        while(p->right)
            p=p->right;
    }
    auto u=new statistic_node<T,V,S>(v,data,p);
    if(v<=p->v)
        p->left=u;
    else p->right=u;
    p->update();
    return rebalance(p);
}

template<typename T, typename S>
statistic_node<T, std::monostate, S>* insert(statistic_node<T, std::monostate, S>* tree,
    const typename std::common_type<T>::type& v,bool or_assign=false)
{
    return insert(tree, v, {},or_assign);
}


/*
* Insert (v,data) into the ordered statistic tree if it does not have a key v
* Otherwise, change the value mapped by v to data.
*/
template<typename T, typename V, typename S>
statistic_node<T, V, S>* insert_or_assign(statistic_node<T, V, S>* tree,const
    typename std::common_type<T>::type& v,const typename std::common_type<V>::type& data)
{
    return insert(tree, v, data, true);
}

/*
* Insert (v,data) into the ordered statistic tree if it does not have a key v
* Otherwise, Do nothing
*/
template<typename T, typename S>
statistic_node<T, std::monostate, S>* insert_or_assign(statistic_node<T, std::monostate, S>* tree, const
    typename std::common_type<T>::type& v)
{
    return insert_or_assign(tree, v, {});
}

/*
* Extract a node from an ordered statistic tree given its key
* The node is not deleted
* @Requirements
* The tree does have a node with the given key
*/
template<typename T,typename V,typename S>
std::pair<statistic_node<T,V,S>*,statistic_node<T,V,S>*> extract(statistic_node<T,V,S>* tree,const 
    typename std::common_type<T>::type& v)
{
    auto p=lower_bound(tree,v);
    if(!p)
        return {nullptr,tree};
    if(!p->left)
    {
        auto w=p->parent?p->parent:p->right;
        if(p->parent)
        {
            if(p->parent->left==p) p->parent->left=p->right;
            else p->parent->right=p->right;
        }
        if(p->right) p->right->parent=p->parent;
        p->right=nullptr;
        p->parent=nullptr;
        return {p,rebalance(w)};
    }
    else if(!p->left->right)
    {
        auto w=p->left;
        if(p->parent)
        {
            if(p->parent->left==p) p->parent->left=p->left;
            else p->parent->right=p->left;
        }
        if(p->right) p->right->parent=w;
        w->right=p->right;
        w->parent=p->parent;
        p->right=nullptr;
        p->left=nullptr;
        p->parent=nullptr;
        return {p,rebalance(w)};
    }
    else
    {
        auto u=p->left;//Position of replacement
        while(u->right)
            u=u->right;
        auto s=u->parent;//Position of path to be updated
        s->right=u->left;
        if(u->left) u->left->parent=s;
        std::swap(u->v,p->v);
	    std::swap(u->data,p->data);
        u->left=nullptr;
        u->right=nullptr;
        u->parent=nullptr;
        return {u,rebalance(s)};
    }

}


template<typename T,typename V,typename S>
statistic_node<T,V,S>* erase(statistic_node<T,V,S>* tree,
    const typename std::common_type<T>::type& v)
{
    auto P=extract(tree,v);
    delete P.first;
    return P.second;
}

template<typename T,typename V, typename S>
statistic_node<T,V,S>* update(statistic_node<T,V,S>*tree, 
    const typename std::common_type<T>::type&  v,
    const typename std::common_type<V>::type& data)
{
    auto p=lower_bound(tree,v);
    p->data=data;
    return rebalance(p);
}

template<typename T,typename V,typename S>
void destroy(statistic_node<T,V,S>*node)
{
    if(!node)
        return;
    destroy(node->left);
    destroy(node->right);
    delete node;
}

/*
* Merge a node with two trees
* @Requirements
* 1. all keys of left are less or equal to the key of the root node
* 2. all keys of right are greater or equal to the key of the root node
* 3. the left & right trees do not have a parent
*/
template<typename T,typename V,typename S>
statistic_node<T, V, S>* merge_with_root(statistic_node<T, V, S>* root, 
    statistic_node<T, V, S> *left, statistic_node<T, V, S> *right)
{
    if (left == nullptr)
        return insert(right, root->v, root->data);
    else if (right == nullptr)
        return insert(left, root->v, root->data);
    auto potential = height(left) - height(right);
    if (potential <=0)
    {
        while (potential < -1)
            potential = height(left) - height(right=right->left);
        if (right && right->parent)
            right->parent->left = root;
        if(right)
            root->parent = right->parent;
    }
    else
    {
        while (potential > 1)
            potential = height(left=left->right) - height(right);
        if (left && left->parent)
            left->parent->right = root;
        if(left)
            root->parent = left->parent;
    }
    root->left = left;
    root->right = right;
    root->update();
    if (left)
        left->parent = root;
    if (right)
        right->parent = root;
    return rebalance(root);
}

/*
* Merge two trees
* @Requirements
* The biggest key of left is smaller than the smallest key of right
*/
template<typename T, typename V, typename S>
statistic_node<T, V, S>* merge(statistic_node<T, V, S>* left, statistic_node<T, V, S>* right)
{
    if (!left)
        return right;
    statistic_node<T, V, S>* last = left;
    while (last->right)
        last = last->right;
    auto [root,L] = extract(last, last->v);
    return merge_with_root(root, L, right);
}
/*
* Split a tree into two trees with respect to threshold
* - The left part is the resultant tree whose keys are smaller than threshold
* - The right part is the resultant tree whose keys are not smaller than threshold
*/
template<typename T, typename V, typename S>
std::pair<statistic_node<T, V, S>*,statistic_node<T,V,S>*> 
    split(statistic_node<T, V, S>* node,T threshold)
{
    statistic_node<T, V, S>* left = nullptr, * right = nullptr;
    if (!node)
        return std::make_pair(left, right);
    if (node->right)
        node->right->parent = nullptr;
    if (node->left)
        node->left->parent = nullptr;
    if (node->v < threshold)
    {
        auto [L, R] = split(node->right, threshold);
        if (L)
            L->parent = nullptr;
        if (R)
            R->parent = nullptr;
        left = merge_with_root(node,node->left,L);
        right = R;
    }
    else
    {
        auto [L, R] = split(node->left, threshold);
        if (L)
            L->parent = nullptr;
        if (R)
            R->parent = nullptr;
        right = merge_with_root(node,R, node->right);
        left = L;
    }
    return std::make_pair(left, right);
}

/*
* Order Statistic:
* It is a Self Balanced Binary Search Tree augmented with an order:
* - The order of an element can be calculated in O(log(n))
* - An element can be selected given its order in O(log(n))
*/

struct order_stats
{
    int size;
    order_stats() {}
    template<typename T,typename V>
    order_stats(T v,V data):size(1){}
    template<typename T,typename V>
    static void update(statistic_node<T,V,order_stats>*node);
};

template<typename T, typename V>
void order_stats::update(statistic_node<T, V, order_stats> *node) {
    node->statistic.size=(node->left?node->left->statistic.size:0)+1+(node->right?node->right->statistic.size:0);
}

template<typename T, typename V,typename OrderStats>
int size(statistic_node<T, V, OrderStats> *node)
{
    return node?node->statistic.size:0;
}

template<typename T,typename V,typename OrderStats>
T order_inf(statistic_node<T,V, OrderStats> *tree,const typename std::common_type<T>::type& v)
{
    if(!tree)
        return -1;
    if(v<tree->v)
        return order_inf(tree->left,v);
    else if(tree->v==v)
    {
        auto o=order_inf(tree->left,v);
        if(o!=-1)
            return o;
        else return size(tree->left);
    }
    else
    {
        auto o=order_inf(tree->right,v);
        if(o!=-1)
            return size(tree->left)+1+o;
        else return -1;
    }
}

template<typename T,typename V,typename OrderStats>
T order_sup(statistic_node<T,V, OrderStats> *tree,const typename std::common_type<T>::type& v)
{
    if(!tree)
        return 0;
    if(v<tree->v)
        return order_sup(tree->left,v);
    else if(tree->v==v)
    {
        if(tree->right && tree->right->v==v)
            return size(tree->left)+1+order_sup(tree->right,v);
        else return size(tree->left);
    }
    else return size(tree->left)+1+order_sup(tree->right,v);
}

template<typename T, typename V, typename OrderStats>
int order(statistic_node<T, V, OrderStats>* tree, const typename std::common_type<T>::type& v )
{
    if (!tree)
        return 0;
    if (v < tree->v)
        return order(tree->left, v);
    else if (tree->v == v)
    {
        if (tree->right && tree->right->v == v)
            return size(tree->left) + 1 + order(tree->right, v);
        else return size(tree->left);
    }
    else return size(tree->left) + 1 + order(tree->right, v);
}

template<typename T,typename V,typename OrderStats>
statistic_node<T,V,OrderStats>* select(statistic_node<T,V, OrderStats> *tree,int o)
{
    int s= size(tree->left);
    if(s==o)
        return tree;
    else if(s<o)
        return select(tree->right,o-s-1);
    else return select(tree->left,o);
}

template<typename T, typename V, typename OrderStats>
statistic_node<T, V, OrderStats>* median(statistic_node<T, V, OrderStats>* tree)
{
    return select(tree, size(tree) / 2);
}

/**
* Sum Statistic:
* It is an Ordered Statistic Tree augmented with a sum acting on data:
* The sum is defined over an associative binary operation having a neutral element
* It supports range sum (L,R) for keys belonging to an interval [L,R[ 
* @Requirements
* 1. a class V
* 2. an associative binary operation O: VxV->V 
* 3. O has a neutral element named 'neutral'. and it is defined as a public static attribute.
* 4. O has an overload for the ternary case: VxVxV->V, and its definition is compatible with the binary case.
* @Notes
* Formally, (V,O) is a monoid
*/
template<typename V, typename O>
struct sum_stats
{
    inline static O F = O();
    int size;
    V sum;
    sum_stats() {}
    template<typename T>
    sum_stats(T v, V data) :size(1), sum(data) {}
    template<typename T>
    static void update(statistic_node<T, V, sum_stats>* node);
    inline static const V& neutral = O::neutral;
};

template<typename V, typename O>
template<typename T>
void sum_stats<V, O>::update(statistic_node<T, V, sum_stats<V, O>>* node) {
    node->statistic.size = (node->left ? node->left->statistic.size : 0) + 1 + (node->right ? node->right->statistic.size : 0);
    node->statistic.sum = F(tree_sum(node->left),node->data, tree_sum(node->right));
}

template<typename T, typename V, typename SumStats>
V tree_sum(statistic_node<T, V, SumStats>* node)
{
    return node ? node->statistic.sum : SumStats::neutral;
}

template<typename T, typename V, typename SumStats>
V prefix_sum(statistic_node<T, V, SumStats>* tree,const  
    typename std::common_type<T>::type& U)
{
    if (!tree)
        return SumStats::neutral;
    else if (tree->v >= U)
        return prefix_sum(tree->left, U);
    else return SumStats::F(tree_sum(tree->left), tree->data, prefix_sum(tree->right, U));
}

template<typename T, typename V, typename SumStats>
V suffix_sum(statistic_node<T, V, SumStats>* tree, const 
    typename std::common_type<T>::type& L)
{
    if (!tree)
        return SumStats::neutral;
    else if (tree->v < L)
        return suffix_sum(tree->right, L);
    else return SumStats::F(suffix_sum(tree->left, L), tree->data, tree_sum(tree->right));
}

template<typename T, typename V, typename SumStats>
V sum(statistic_node<T, V, SumStats>* tree, const 
    typename std::common_type<T>::type& L, const typename std::common_type<T>::type& R)
{
    if (!tree)
        return SumStats::neutral;
    if (tree->v < L)
        return sum(tree->right, L, R);
    else if (tree->v >= R)
        return sum(tree->left, L, R);
    else return SumStats::F(suffix_sum(tree->left, L), tree->data, prefix_sum(tree->right, R));
}

template<typename T, typename V, typename SumStats>
V prefix_index_sum(statistic_node<T, V, SumStats>* tree, int n)
{
    if (!tree || n <= 0)
        return SumStats::neutral;
    else if (size(tree->left) >= n)
        return prefix_index_sum(tree->left, n);
    else return SumStats::F(tree_sum(tree->left), tree->data, prefix_index_sum(tree->right, n - 1 - size(tree->left)));
}

template<typename T, typename V, typename SumStats>
V suffix_index_sum(statistic_node<T, V, SumStats>* tree, int n)
{
    if (!tree)
        return SumStats::neutral;
    else if (size(tree->right) >= n)
        return suffix_index_sum(tree->right, n);
    else return SumStats::F(suffix_index_sum(tree->left, n - 1 - size(tree->right)), tree->data, tree_sum(tree->right));
}

template<typename T, typename V, typename SumStats>
V index_sum(statistic_node<T, V, SumStats>* tree, int a, int b)
{
    if (!tree || a >= b)
        return SumStats::neutral;
    int left_size = size(tree->left);
    if (left_size < a)
        return SumStats::F(tree->data,index_sum(tree->right, a - left_size - 1, b - left_size - 1));
    else if (left_size >= b)
        return index_sum(tree->left, a, b);
    else return SumStats::F(suffix_index_sum(tree->left, left_size), tree->data, prefix_index_sum(tree->right, left_size - 1));
}

/** Key Sum Statistic:
* It is an Ordered Statistic Tree augmented with a sum acting on keys:
* The sum is defined over an associative binary operation having a neutral element
* It supports range sum (L,R) for keys belonging to the interval [L,R[ 
* @Requirements
* 1. a class T
* 2. an associative binary operation O: TxT->T 
* 3. O has a neutral element named 'neutral'. and it is defined as a public static attribute.
* 4. O has an overload for the ternary case: TxTxT->T, and its definition is compatible with the binary case.
* @Notes
* Formally, (T,O) is a monoid
* The order of T does not need to be compatible with O
*/

template<typename T, typename O>
struct key_sum_stats
{
    inline static O F = O();
    int size;
    T key_sum;
    key_sum_stats() {}
    template<typename V>
    key_sum_stats(T v, V data) :size(1), key_sum(v) {}
    template<typename V>
    static void update(statistic_node<T, V, key_sum_stats>* node);
    inline static const T& key_neutral = O::neutral;
};

template<typename T, typename O>
template<typename V>
void key_sum_stats<T, O>::update(statistic_node<T, V, key_sum_stats<T, O>>* node) {
    node->statistic.size = (node->left ? node->left->statistic.size : 0) + 1 + (node->right ? node->right->statistic.size : 0);
    node->statistic.key_sum = F(tree_key_sum(node->left), node->v, tree_key_sum(node->right));
}


template<typename T, typename V, typename KeySumStats>
T tree_key_sum(statistic_node<T, V, KeySumStats>* node)
{
    return node ? node->statistic.key_sum : KeySumStats::key_neutral;
}

template<typename T, typename V, typename KeySumStats>
T prefix_key_sum(statistic_node<T, V, KeySumStats>* tree, const 
    typename std::common_type<T>::type &U)
{
    if (!tree)
        return KeySumStats::key_neutral;
    else if (tree->v >= U)
        return prefix_key_sum(tree->left, U);
    else return KeySumStats::F(tree_key_sum(tree->left), tree->v, prefix_key_sum(tree->right, U));
}

template<typename T, typename V, typename KeySumStats>
T suffix_key_sum(statistic_node<T, V, KeySumStats>* tree, const 
    typename std::common_type<T>::type &L)
{
    if (!tree)
        return KeySumStats::key_neutral;
    else if (tree->v < L)
        return suffix_key_sum(tree->right, L);
    else return KeySumStats::F(suffix_key_sum(tree->left, L), tree->v, tree_key_sum(tree->right));
}

template<typename T, typename V, typename KeySumStats>
T key_sum(statistic_node<T, V, KeySumStats>* tree, const typename std::common_type<T>::type& L,
    const typename std::common_type<T>::type& R)
{
    if (!tree)
        return KeySumStats::key_neutral;
    if (tree->v < L)
        return key_sum(tree->right, L, R);
    else if (tree->v >= R)
        return key_sum(tree->left, L, R);
    else return KeySumStats::F(suffix_key_sum(tree->left, L), tree->v, prefix_key_sum(tree->right, R));
}

template<typename T, typename V, typename KeySumStats>
T prefix_index_key_sum(statistic_node<T, V, KeySumStats>* tree, int n)
{
    if (!tree || n <= 0)
        return KeySumStats::key_neutral;
    else if (size(tree->left) >= n)
        return prefix_index_key_sum(tree->left, n);
    else return KeySumStats::F(tree_key_sum(tree->left), tree->v, prefix_index_key_sum(tree->right, n - 1 - size(tree->left)));
}

template<typename T, typename V, typename KeySumStats>
T suffix_index_key_sum(statistic_node<T, V, KeySumStats>* tree, int n)
{
    if (!tree)
        return KeySumStats::key_neutral;
    else if (size(tree->right) >= n)
        return suffix_index_key_sum(tree->right, n);
    else return KeySumStats::F(suffix_index_key_sum(tree->left, n-1-size(tree->right)), tree->v, tree_key_sum(tree->right));
}

template<typename T, typename V, typename KeySumStats>
T index_key_sum(statistic_node<T, V, KeySumStats>* tree, int a,int b)
{
    if (!tree || a>=b)
        return KeySumStats::key_neutral;
    int left_size = size(tree->left);
    if (left_size < a)
        return KeySumStats::F(tree->v,index_key_sum(tree->right, a-left_size-1, b-left_size-1));
    else if (left_size >= b)
        return index_key_sum(tree->left, a, b);
    else return KeySumStats::F(suffix_index_key_sum(tree->left, left_size-a), tree->v, prefix_index_key_sum(tree->right, b-left_size-1));
}

/*
* The following aliases gives the possibility to a much shorter and compact code
* For example:
* 1. sum_node_t<int,double,multiplies_t> is a shorthand for this atrocity:
*   statistic_node<int,double,sum_stats<double,multiplies_t<double>>>
*   this is a statistic_node whose key is of type int and values are doubles, the considered statistic is 
*   an aggregation stats whose operation is multiplication
* The first two considers the case when the binary operator's class is a template
* The last two considers the case when the binary operator's class is not a template
*/

//
template<typename T, typename V, template<typename S = T> typename O>
using sum_node_t = statistic_node<T, V, sum_stats<V, O<V>>>;

template<typename T, template<typename S = T> typename O, typename V=std::monostate>
using key_sum_node_t = statistic_node<T, V, key_sum_stats<T, O<T>>>;

template<typename T, typename V, typename O>
using sum_node = statistic_node<T, V, sum_stats<V, O>>;

template<typename T, typename O, typename V = std::monostate>
using key_sum_node = statistic_node<T, V, key_sum_stats<T, O>>;

template<typename T, typename V = std::monostate>
using order_node = statistic_node<T, V, order_stats>;

/*
* A
*/
#endif

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

template<typename O>
struct prefix_array
{
    using R=typename O::type;
    using type=typename O::type;
    std::vector<R> A;
    std::vector<R> P;
    inline static constexpr O F=O();
    prefix_array(const std::vector<R> &_A):A(_A),P(_A.size()+1)
    {
        P[0]=O::neutral;
        for(int i=0;i<A.size();i++)
            P[i+1]=F(P[i],A[i]);
    }

    R query(int l,int r)
    {
        return F(F.inv(P[l]),P[r]);
    }

    void update(int i,R u)
    {
        A[i]=u;
        for(int j=i+1;j<P.size();j++)
            P[j]=F(P[j-1],A[j-1]);
    }
};

template<typename O>
struct segment_tree
{
    using R=typename O::type;
    using type=R;
    std::vector<std::vector<R>> S;
    std::vector<R> A;
    int n,h;
    segment_tree(const std::vector<R> &_A):A(_A)
    {
        n=bit_ceil(A.size());
        A.resize(n,O::neutral);
        int m=n;
        h=0;
        while(m)
        {
            m/=2;
            h++;
        }
        S.resize(h);
        for(int i=0;i<h;i++)
            S[i].resize(1<<i);
        build();
    }

    void update(int i,R u)
    {
        A[i]=u;
        S[h-1][i]=u;
        int m=h-2;
        i/=2;
        while(m>=0)
        {
            S[m][i]=F(S[m+1][2*i],S[m+1][2*i+1]);
            m--;
            i/=2;
        }
    }

    R query(int l,int r)
    {
        return query(std::max(l,0),std::min(r,n),0,n,0);
    }
private:
    inline static O F=O();
    void build()
    {
        for(int i=0;i<n;i++)
            S.back()[i]=A[i];
        for(int i=h-2;i>=0;i--) for(int k=0;k<(1<<i);k++)
                S[i][k]=F(S[i+1][2*k],S[i+1][2*k+1]);
    }
    R query(int l,int r,int a,int b,int depth)
    {
        if(l>=r)
            return O::neutral;
        if(l==a && r==b)
            return S[depth][l>>(h-1-depth)];
        int mid=(a+b)/2;
        if(mid>r)
            return query(l,r,a,mid,depth+1);
        else if(mid<l)
            return query(l,r,mid,b,depth+1);
        else
            return F(query(l,mid,a,mid,depth+1),query(mid,r,mid,b,depth+1));
    }
};

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

template<typename O>
struct segment_matrix
{
    using R=typename O::type;
    using type=R;
    std::vector<std::vector<segment_tree<O>>> S;
    std::vector<segment_tree<O>> segment_forest;
    std::vector<std::vector<R>> A;
    int n,h;
    segment_matrix(std::vector<std::vector<R>> &&_A):A(std::move(_A)),h(0)
    {
        n=bit_ceil(A.size());
        A.resize(n,std::vector<R>(bit_ceil(A[0].size()),O::neutral));
        int w=n;
        while(w)
        {
            w/=2;
            h++;
        }
        S.resize(h);
        int m=A[0].size();
        for(int i=0;i<m;i++)
        {
            std::vector<R> C;
            for(int j=0;j<n;j++)
                C.push_back(A[j][i]);
            segment_forest.emplace_back(C);
        }
        std::vector<R> C(m);
        for(int i=h-1;i>=0;i--) for(int j=0;j<(1<<i);j++)
            {
                for(int p=0;p<m;p++)
                    C[p]=segment_forest[p].S[i][j];
                S[i].emplace_back(C);
            }
    }

    void update(int i,int j,R u)
    {
        A[i][j]=u;
        segment_forest[j].update(i,u);
        int r=h-1;
        while(r>=0)
        {
            S[r][i].update(j,segment_forest[j].S[r][i]);
            r--;
            i/=2;
        }
    }

    R query(int l,int r,int p,int q)
    {
        return query(l,r,p,q,0,n,0);
    }
private:
    inline static O F=O();

    R query(int l,int r,int p,int q,int a,int b,int depth)
    {
        if(l>=r)
            return O::neutral;
        if(l==a && r==b)
            return S[depth][l>>(h-1-depth)].query(p,q);
        int mid=(a+b)/2;
        if(mid>r)
            return query(l,r,p,q,a,mid,depth+1);
        else if(mid<l)
            return query(l,r,p,q,mid,b,depth+1);
        else
            return F(query(l,mid,p,q,a,mid,depth+1),query(mid,r,p,q,mid,b,depth+1));
    }

};

template<typename T,typename O>
struct fenwick_tree {
    int n;
    std::vector<T> bit;
    inline static O F = O();

    fenwick_tree(int _n):n(_n),bit(n,O::neutral){}
    fenwick_tree(const std::vector<T> &X) : fenwick_tree(X.size())
    {
        for(int i=0;i<n;i++)
            update(i,X[i]);
    }
    T sum(int x) {
        if(x<0)
            return O::neutral;
        T ret = O::neutral;
        for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
            ret = F(ret,bit[i]);
        return ret;
    }

    T query(int a,int b)
    {
        return F(F.inv(sum(a-1)),sum(b));
    }

    T sum(int a,int b)
    {
        return query(a,b);
    }

    void add(int x, T delta) {
        for (int i = x; i < n; i = i | (i + 1))
            bit[i] = F(bit[i], delta);
    }

    void update(int x, T delta) {
        add(x,F(F.inv(sum(x,x)),delta));
    }
};


template<typename T,typename O>
struct fenwick_matrix {
    inline static O F=O();
    int n, m;
    std::vector<std::vector<T>> bit;

    fenwick_matrix(int _n,int _m):n(_n),m(_m),bit(n,std::vector<T>(m,O::neutral)){}
    int sum(int x, int y) {
        if(x<0||y<0)
            return O::neutral;
        int ret = 0;
        for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
            for (int j = y; j >= 0; j = (j & (j + 1)) - 1)
                ret = F(ret,bit[i][j]);
        return ret;
    }

    int sum(int a,int b,int c,int d)
    {
        //To Do
        //........................................
    }


    int query(int a,int b,int c,int d)
    {
        return F(F.inv(sum(a,c),sum(b,d)));
    }

    void add(int x, int y, int delta) {
        for (int i = x; i < n; i = i | (i + 1))
            for (int j = y; j < m; j = j | (j + 1))
                bit[i][j] = F(bit[i][j],delta);
    }

    void update(int x, int y, int delta) {
        add(x,y,F(F.inv(sum(x,x,y,y)),delta));
    }
};


template<typename T,typename O>
struct sparse_segment_tree
{
    sum_node<int, T, O>* tree;
public:
    sparse_segment_tree():tree(nullptr) {}
    ~sparse_segment_tree()
    {
        destroy(tree);
        tree = nullptr;
    }

    void insert(int k, const T& v)
    {
        tree = ::insert(tree, k, v);
    }

    void update(int k, const T& v)
    {
        tree = insert_or_assign(k, v);
    }

    void erase(int k)
    {
        tree = ::erase(tree, k);
    }

    T query(int l, int r) const
    {
        return sum(tree, l, r);
    }

    T index_query(int l, int r) const
    {
        return index_sum(tree,l, r);
    }
};

template<typename T, typename O>
class dynamic_segment_tree
{
    sum_node<int, T, O>* tree;
    int size;
public:
    dynamic_segment_tree() :tree(nullptr),size(0) {}
    ~dynamic_segment_tree()
    {
        destroy(tree);
        tree = nullptr;
        size = 0;
    }

    void push_back(const T& v)
    {
        tree = insert(tree, size++, v);
    }

    void update(int k, const T& v)
    {
        tree = insert_or_assign(k, v);
    }

    void pop_back()
    {
        tree = erase(tree, size--);
    }

    T query(int l, int r) const
    {
        return sum(tree, l, r);
    }

    T index_query(int l, int r) const
    {
        return query(l, r);
    }
};

template<typename T, typename O>
struct ordered_segment_tree
{
    key_sum_node<T,O>* tree;
public:
    ordered_segment_tree() :tree(nullptr) {}
    ~ordered_segment_tree()
    {
        destroy(tree);
        tree = nullptr;
    }

    void insert(const T& v)
    {
        tree = ::insert(tree, v);
    }

    void erase(const T&v)
    {
        tree = ::erase(tree, v);
    }

    T query(const T& l, const T& r) const
    {
        return key_sum(tree, l, r);
    }

    T index_query(int a, int b) const
    {
        return index_key_sum(tree, a, b);
    }
};


#endif
#include <memory>
#include <queue>
#include <random>
#include <chrono>

/*
 * EdgeIdentifier is a struct that identifies an edge in the tree
 * by its two endpoints and its index in the adjacency list of the first endpoint
 * */
struct EdgeIdentifier
{
    int u,v;
    std::strong_ordering operator<=>(const EdgeIdentifier &rhs) const noexcept
    {
        return std::tie(u,v)<=>std::tie(rhs.u,rhs.v);
    }
    bool operator==(const EdgeIdentifier &rhs) const noexcept
    {
        return std::tie(u,v)==std::tie(rhs.u,rhs.v);
    }
    EdgeIdentifier(int _u,int _v):u(_u),v(_v){}
};

template<>
class std::hash<EdgeIdentifier>
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    inline static constexpr integer M=1e9+7;
    std::uniform_int_distribution<integer> dis{1,M-1};
    integer x=dis(gen),y=dis(gen),z=dis(gen);
public:
    std::size_t operator()(const EdgeIdentifier &Z) const noexcept
    {
        return (Z.u*x+Z.v*y+z)%M;
    }
};

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
        lca_data=std::make_unique<sparse_array<min_t<HeightData>>>(A);
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


protected:
    using HeightData=std::pair<int,int>;
    using EnpointsData = std::pair<int,int>;
    std::unique_ptr<sparse_array<min_t<HeightData>>> lca_data;
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
        lca_data=std::make_unique<sparse_array<min_t<HeightData>>>(A);
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
    std::unique_ptr<sparse_array<min_t<HeightData>>> lca_data;
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


#endif //CPLIBRARY_TREE_H


template<typename O, typename RMQ>
struct HeavyLightTree : public WeightedTree<typename O::type>
{
    using type=typename O::type;
    using Weight=type;
private:
public:
    using WeightedTree<Weight>::HLD;
    using WeightedTree<Weight>::WeightedTree;
    using AdjacentType=WeightedTree<Weight>::AdjacentType;
    using WeightedTree<Weight>::root;
    using WeightedTree<Weight>::parent;
    inline static O F{};

    void buildStatistics(TreeStats stats=TreeStats::RANGE_QUERIES)
    {
        WeightedTree<Weight>::buildStatistics(stats);
        if(stats & TreeStats::RANGE_QUERIES) buildRangeStatistics();
    }

    Weight query(int u,int v)
    {
        auto lca=WeightedTree<Weight>::leastCommonAncestor(u,v);
        return F(query_with_lca(u,lca,false),query_with_lca(v,lca,true));
    }

    void update(int u,int v, Weight w)
    {
        if(parent[u] && parent[u]->first==v)
            return update(u,w);
        if(parent[v] && parent[v]->first==u)
            return update(v,w);
        throw std::invalid_argument("u and v must be adjacent");
    }

    void update(int u,Weight w)
    {
        if(!parent[u])
            throw std::invalid_argument("u must not be the root");
        parent[u]->second=w;
        auto v=parent[u]->first;
        if(WeightedTree<Weight>::HLD.is_heavy[u])
        {
            auto &R_left=S_left[HLD.HLD_mapper[v].hld_id];
            auto &R_right=S_right[HLD.HLD_mapper[v].hld_id];
            R_left.update(HLD.component_size[HLD.HLD_mapper[v].hld_id] - HLD.HLD_mapper[v].index - 1, w);
            R_right.update(HLD.HLD_mapper[v].index, w);
        }
    }

protected:


    //TODO: Support for non-commutative operations
    Weight query_with_lca(int u,int lca, bool invert)
    {
        int o=u;
        if(u==lca)
            return O::neutral;
        Weight R=O::neutral;
        while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id)
        {
            auto [x,y]=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id];
            auto b=HLD.HLD_mapper[u].index;
            auto a= HLD.HLD_mapper[x].index;
            auto &S=S_right;
            if(!invert)
                R=F(R,S[HLD.HLD_mapper[u].hld_id].query(a, b));
            else
                R=F(S[HLD.HLD_mapper[u].hld_id].query(a, b), R);
            u=x;
            while(u!=lca && !WeightedTree<Weight>::HLD.is_heavy[u])
            {
                if(!invert) R=F(R,parent[u]->second);
                else R=F(parent[u]->second,R);
                u=parent[u]->first;
            }
        }
        auto b=HLD.HLD_mapper[u].index;
        auto a= HLD.HLD_mapper[lca].index;
        auto &S=S_right;
        if(!invert)
            R=F(R,S[HLD.HLD_mapper[u].hld_id].query(a, b));
        else
            R=F(S[HLD.HLD_mapper[u].hld_id].query(a, b), R);
        return R;
    }

    std::vector<RMQ> S_left, S_right;

    void buildRangeStatistics()
    {
        for(auto &C:HLD.components)
        {
            HLD.component_size.push_back(C.size());
            if(C.empty()) C.emplace_back(O::neutral);
            S_right.emplace_back(C);
            std::reverse(C.begin(),C.end());
            S_left.emplace_back(std::move(C));
        }
    }
};

template<typename O, typename RMQ>
struct HeavyLightNodeTree : protected HeavyLightTree<O,RMQ>
{
    using type=typename O::type;
    using Weight=type;
protected:
    void buildRepresentation(int root)
    {
        reRoot(root);
        for(int u=0;u<n;u++) for(auto &[v,w]: this->reverseList[u])
            {
                w=node_weight[v];
                parent[v]=std::make_pair(u,w);
                this->adjacencyList[v].front().second=w;
            }
    }
    using WeightedTree<Weight>::children;
    using HeavyLightTree<O,RMQ>::query_with_lca;
    using HeavyLightTree<O,RMQ>::F;
    using WeightedTree<Weight>::reRoot;
public:
    std::vector<Weight> node_weight;
    HeavyLightNodeTree(int n,int root):HeavyLightTree<O,RMQ>(n,root),node_weight(n,O::neutral){}
    HeavyLightNodeTree(int n):HeavyLightNodeTree(n,0){}
    using AdjacentType=WeightedTree<Weight>::AdjacentType;
    using WeightedTree<Weight>::root;
    using WeightedTree<Weight>::parent;
    using WeightedTree<Weight>::n;

    void setWeight(int u,Weight w)
    {
        node_weight[u]=w;
    }

    void connect(int u,int v)
    {
        HeavyLightTree<O,RMQ>::connect(u,v,O::neutral);
    }

    void setParent(int u,int v)
    {
        HeavyLightTree<O,RMQ>::setParent(u,v,O::neutral);
    }

    void build(int root=0,TreeStats stats= TreeStats::RANGE_QUERIES)
    {
        buildRepresentation(root);
        HeavyLightTree<O,RMQ>::buildStatistics(stats);
    }

    Weight query(int u,int v)
    {
        auto lca=WeightedTree<Weight>::leastCommonAncestor(u,v);
        return F(query_with_lca(u,lca,false),node_weight[lca],query_with_lca(v,lca,true));
    }

    void update(int u, Weight w)
    {
        node_weight[u]=w;
        if(parent[u])
            HeavyLightTree<O, RMQ>::update(u, w);
    }

    using HeavyLightTree<O,RMQ>::distance;
};


template<typename O, typename RMQ>
struct CommutativeHeavyLightTree : public WeightedTree<typename O::type>
{
    using Weight=typename O::type;
    using type=Weight;
private:
public:
    using WeightedTree<Weight>::HLD;
    using WeightedTree<Weight>::WeightedTree;
    using AdjacentType=WeightedTree<Weight>::AdjacentType;
    using WeightedTree<Weight>::root;
    using WeightedTree<Weight>::parent;
    inline static O F{};

    void buildStatistics(TreeStats stats=TreeStats::RANGE_QUERIES)
    {
        WeightedTree<Weight>::buildStatistics(stats);
        if(stats & TreeStats::RANGE_QUERIES) buildRangeStatistics();
    }

    Weight query(int u,int v)
    {
        auto lca=WeightedTree<Weight>::leastCommonAncestor(u,v);
        return F(query_with_lca(u,lca),query_with_lca(v,lca));
    }

    void update(int u,int v, Weight w)
    {
        if(parent[u] && parent[u]->first==v)
            return update(u,w);
        if(parent[v] && parent[v]->first==u)
            return update(v,w);
        throw std::invalid_argument("u and v must be adjacent");
    }

    void update(int u,Weight w)
    {
        if(!parent[u])
            throw std::invalid_argument("u must not be the root");
        parent[u]->second=w;
        auto v=parent[u]->first;
        if(WeightedTree<Weight>::HLD.is_heavy[u])
        {
            auto &R=S[HLD.HLD_mapper[v].hld_id];
            R.update(HLD.HLD_mapper[v].index, w);
        }
    }

protected:
    std::unordered_map<EdgeIdentifier,Weight> cache;
    Weight query_with_lca(int u,int lca)
    {
        if(u==lca)
            return O::neutral;
        Weight R=O::neutral;
        while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id)
        {
            auto [x,y]=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id];
            auto b=HLD.HLD_mapper[u].index;
            auto a= HLD.HLD_mapper[x].index;
            R=F(S[HLD.HLD_mapper[u].hld_id].query(a, b), R);
            u=x;
            while(u!=lca && !WeightedTree<Weight>::HLD.is_heavy[u])
            {
                R=F(R,parent[u]->second);
                u=parent[u]->first;
            }
        }
        auto b=HLD.HLD_mapper[u].index;
        auto a= HLD.HLD_mapper[lca].index;
        R=F(R,S[HLD.HLD_mapper[u].hld_id].query(a, b));
        return R;
    }

    std::vector<RMQ> S;

    void buildRangeStatistics()
    {
        for(auto &C:HLD.components)
        {
            HLD.component_size.push_back(C.size());
            if(C.empty()) C.emplace_back(O::neutral);
            S.emplace_back(C);
        }
    }
};



template<typename O, typename RMQ>
struct CommutativeHeavyLightNodeTree : protected CommutativeHeavyLightTree<O,RMQ>
{
    using type=typename O::type;
    using Weight=type;
protected:
    void buildRepresentation(int root)
    {
        reRoot(root);
        for(int u=0;u<n;u++) for(auto &[v,w]: this->reverseList[u])
            {
                w=node_weight[v];
                parent[v]=std::make_pair(u,w);
                this->adjacencyList[v].front().second=w;
            }
    }
    using WeightedTree<Weight>::children;
    using CommutativeHeavyLightTree<O,RMQ>::query_with_lca;
    using CommutativeHeavyLightTree<O,RMQ>::F;
    using WeightedTree<Weight>::reRoot;
public:
    std::vector<Weight> node_weight;
    CommutativeHeavyLightNodeTree(int n,int root): CommutativeHeavyLightTree<O,RMQ>(n,root),node_weight(n,O::neutral){}
    CommutativeHeavyLightNodeTree(int n):CommutativeHeavyLightNodeTree(n,0){}
    using AdjacentType=WeightedTree<Weight>::AdjacentType;
    using WeightedTree<Weight>::root;
    using WeightedTree<Weight>::parent;
    using WeightedTree<Weight>::n;

    void setWeight(int u,Weight w)
    {
        node_weight[u]=w;
    }

    void connect(int u,int v)
    {
        CommutativeHeavyLightTree<O,RMQ>::connect(u,v,O::neutral);
    }

    void setParent(int u,int v)
    {
        CommutativeHeavyLightTree<O,RMQ>::setParent(u,v,O::neutral);
    }

    void build(int root=0,TreeStats stats= TreeStats::RANGE_QUERIES)
    {
        buildRepresentation(root);
        CommutativeHeavyLightTree<O,RMQ>::buildStatistics(stats);
    }

    Weight query(int u,int v)
    {
        auto lca=WeightedTree<Weight>::leastCommonAncestor(u,v);
        return F(query_with_lca(u,lca),node_weight[lca],query_with_lca(v,lca));
    }

    void update(int u, Weight w)
    {
        node_weight[u]=w;
        if(parent[u])
            CommutativeHeavyLightTree<O, RMQ>::update(u, w);
    }

    using CommutativeHeavyLightTree<O,RMQ>::distance;
};
#endif //CPLIBRARY_RANGE_QUERIES_H
#include <chrono>

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int n,q;
    std::cin >> n >> q;
    CommutativeHeavyLightNodeTree<max_t<integer>,segment_tree<max_t<integer>>> H(n);
    std::vector<integer> A(n);
    for(auto &a:A)
        std::cin >> a;
    for(int i=0;i<n;i++)
        H.setWeight(i,A[i]);
    for(int i=0;i<n-1;i++)
    {
        int u,v;
        std::cin >> u >> v;
        u--;
        v--;
        H.connect(u,v);
    }
    H.build(0);
    auto t3=std::chrono::high_resolution_clock::now();
    for(int i=0;i<q;i++)
    {
        int b;
        std::cin >> b;
        if(b==1)
        {
            int s;
            integer x;
            std::cin >> s >> x;
            s--;
            H.update(s,x);
        }
        else
        {
            int u,v;
            std::cin >> u >> v;
            u--;
            v--;
            std::cout << H.query(u,v) << '\n';
        }
    }
    auto t4=std::chrono::high_resolution_clock::now();
    //std::cerr << "dt: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << std::endl;
}
