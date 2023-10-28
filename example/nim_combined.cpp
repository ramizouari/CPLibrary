#include <stdexcept>
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
    std::variant<std::vector<T>,std::span<T>> value;
    view_or_value(std::vector<T> _value):value(std::move(_value)){}
    view_or_value(std::span<T> _value):value(std::move(_value)){}
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
//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_GRAPH_ALGORITHMS_H
#define CPLIBRARY_GRAPH_ALGORITHMS_H
//
// Created by ramizouari on 27/10/23.
//

#ifndef CPLIBRARY_BELLMAN_FORD_H
#define CPLIBRARY_BELLMAN_FORD_H
#ifndef __ORDER_H__
#define __ORDER_H__
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
#include <compare>
#include <variant>

/*
* Let (S,<=) be a totally ordered set.
* By definition, a closure order of (S,<=) is a totally ordered set (S',<=)
* Where S'=S∪{a,b} where a,b are not members of S
* Furthermore, we define:
* 1. a<=s for all s in S
* 2. s<=b for all s in S
*/
struct inf_minus_t;

struct inf_t : public std::monostate
{
};

constexpr struct inf_plus_t :public inf_t
{
    std::strong_ordering operator<=>(const inf_plus_t&) const = default;
    bool operator==(const inf_plus_t&) const = default;
    bool operator==(const inf_minus_t&) const
    {
        return false;
    }

    std::strong_ordering operator<=>(const inf_minus_t&) const
    {
        return 1 <=> 0;
    }
} inf, inf_plus;

constexpr struct inf_minus_t: public inf_t
{
    bool operator==(const inf_plus_t&) const
    {
        return false;
    }
    std::strong_ordering operator<=>(const inf_plus_t&) const
    {
        return 0 <=> 1;
    }
    std::strong_ordering operator<=>(const inf_minus_t&) const = default;
} inf_min;

constexpr inf_minus_t operator-(const inf_plus_t&)
{
    return inf_min;
}

constexpr inf_plus_t operator-(const inf_minus_t&)
{
    return inf_plus;
}

template<typename S>
using order_closure = std::variant<inf_minus_t, S, inf_plus_t>;
using extended_real = order_closure<real>;
using extended_integer = order_closure<integer>;

/*
* Algebraic Operations on an order closure
* Formally, if S has also a group or ring like structure, we can augment the definition of addition, multiplication
* on almost all elements of S'.
* However, S' does not have the algebraic structure of S
*/
template<typename S>
order_closure<S> operator-(const order_closure<S>& A)
{
    return std::visit([](const auto& B)->order_closure<S> {return -B; }, A);
}

template<typename S>
order_closure<S> operator+(const order_closure<S>& A, const order_closure<S>& B)
{
    if (A.index() == 1 && B.index() != 1)
        return B;
    else if (A.index() != 1 && B.index() == 1)
        return A;
    else if (A.index() == 1 && B.index() == 1)
        return std::get<S>(A) + std::get<S>(B);
    else if (A.index() == B.index())
        return A;
    else return S{};
}

template<typename S>
order_closure<S> operator-(const order_closure<S>& A, const order_closure<S>& B)
{
    if (A.index() == 1 && B.index() != 1)
        return -B;
    else if (A.index() != 1 && B.index() == 1)
        return A;
    else if (A.index() == 1 && B.index() == 1)
        return std::get<S>(A) - std::get<S>(B);
    else if (A.index() == B.index())
        return S{};
    else return A;
}

template<typename S>
order_closure<S> operator-(const order_closure<S>& A, const S& k)
{
    if (A.index() == 1)
        return std::get<S>(A) - k;
    else return A;
}

template<typename S>
order_closure<S> operator-(const S& k,const order_closure<S>& A)
{
    if (A.index() == 1)
        return std::get<S>(A) - k;
    else return -A;
}

template<typename S>
order_closure<S> operator+(const order_closure<S>& A, const S& k)
{
    if (A.index() == 1)
        return std::get<S>(A) + k;
    else return A;
}

template<typename S>
order_closure<S> operator+(const S& k, const order_closure<S>& A)
{
    return A + k;
}

template<typename S>
order_closure<S> operator*(const order_closure<S>& A, const order_closure<S>& B)
{
    if (A.index() == 1 && B.index() != 1)
    {
        if (std::get<S>(A) == 0)
            return 0;
        else if (std::get<S>(A) > 0)
            return B;
        else return -B;
    }
    else if (A.index() != 1 && B.index() == 1)
    {
        if (std::get<S>(B) == 0)
            return 0;
        else if (std::get<S>(B) > 0)
            return A;
        else return -A;
    }
    else if (A.index() == 1 && B.index() == 1)
        return std::get<S>(A) * std::get<S>(B);
    return A.index() == B.index() ? order_closure<S>(inf) : order_closure<S>(-inf);
}

template<typename S>
order_closure<S> operator*(const order_closure<S>& A, const S& k)
{
    if (A.index() == 1)
        return std::get<S>(A) * k;
    else if (k == 0)
        return 0;
    else if (k > 0)
        return A;
    else return -A;
}

template<typename S>
order_closure<S> operator*(const S& k, const order_closure<S>& A)
{
    return operator*(A, k);
}

template<typename S>
order_closure<S> operator/(const order_closure<S>& A, const order_closure<S>& B)
{
    if (A.index() == 1 && B.index() != 1)
        return 0;
    else if (A.index() != 1 && B.index() == 1)
    {
        if (std::get<S>(B) >= 0)
            return A;
        else return -A;
    }
    else if (A.index() == 1 && B.index() == 1)
        return std::get<S>(A) / std::get<S>(B);
    return A.index() == B.index() ? 1 : -1;
}

template<typename S>
order_closure<S> operator/(const order_closure<S>& A, const S& k)
{
    if (A.index() == 1)
        return std::get<S>(A) / k;
    else if (k >= 0)
        return A;
    else return -A;
}

template<typename S>
order_closure<S> operator/(const S& k, const order_closure<S>& A)
{
    if (A.index() != 1)
        return 0;
    else return k / std::get<S>(A);
}

template<typename S>
order_closure<S>& operator+=(order_closure<S>& A, const order_closure<S>& B)
{
    return A = A + B;
}

template<typename S>
order_closure<S>& operator-=(order_closure<S>& A, const order_closure<S>& B)
{
    return A = A - B;
}

template<typename S>
order_closure<S>& operator*=(order_closure<S>& A, const order_closure<S>& B)
{
    return A = A * B;
}

template<typename S>
order_closure<S>& operator/=(order_closure<S>& A, const order_closure<S>& B)
{
    return A = A / B;
}

template<typename S>
order_closure<S>& operator+=(order_closure<S>& A, const S& B)
{
    return A = A + B;
}

template<typename S>
order_closure<S>& operator-=(order_closure<S>& A, const S& B)
{
    return A = A - B;
}

template<typename S>
order_closure<S>& operator*=(order_closure<S>& A, const S& B)
{
    return A = A * B;
}

template<typename S>
order_closure<S>& operator/=(order_closure<S>& A, const S& B)
{
    return A = A / B;
}

template<typename S>
using base_order_type = decltype(std::declval<S>() <=> std::declval<S>());

template<typename S>
base_order_type<S> operator<=>(const order_closure<S>& A, const S&B)
{
    using order_type=base_order_type<S>;
    if (A.index() == 1)
        return std::get<S>(A) <=> B;
    else if (A.index() == 0)
        return order_type::less;
    else return order_type::greater;
}

template<typename S>
bool operator==(const order_closure<S>& A, const S&B)
{
    return A <=> B == 0;
}

template<typename S>
base_order_type<S> operator<=>(const order_closure<S>&A,inf_minus_t B)
{
    using order_type=base_order_type<S>;
    return A.index() == 0 ? order_type::equivalent : order_type::greater;
}

template<typename S>
base_order_type<S> operator<=>(const order_closure<S>&A,inf_plus_t B)
{
    using order_type=base_order_type<S>;
    return A.index() == 2 ? order_type::equivalent : order_type::less;
}

template<typename S>
bool operator==(const order_closure<S>&A,inf_minus_t B)
{
    return A <=> B == 0;
}

template<typename S>
bool operator==(const order_closure<S>&A,inf_plus_t B)
{
    return A <=> B == 0;
}

template<typename S>
std::ostream &operator<<(std::ostream &os, const order_closure<S> &A)
{
    if (A.index() == 1)
        os << std::get<S>(A);
    else if (A.index() == 0)
        os << "-inf";
    else os << "inf";
    return os;
}
#endif
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

/**
 * @brief This the class of directed Graphs
 * @details Each Graph is a couple G=(V,E) where V={0,...,n-1} are nodes, and E is a subset of VxV
 * @param n the size of the Graph
 * */

namespace graph
{
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

}

#endif //CPLIBRARY_GRAPH_H
namespace graph::algorithms
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
//
// Created by ramizouari on 27/10/23.
//

#ifndef CPLIBRARY_DIJKSTRA_H
#define CPLIBRARY_DIJKSTRA_H
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
//
// Created by ramizouari on 27/10/23.
//

#ifndef CPLIBRARY_BFS_H
#define CPLIBRARY_BFS_H
#include <queue>

namespace graph::algorithms
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

#endif //CPLIBRARY_ALGORITHMS_H
#include <iostream>

struct CollatzGraph:public graph::AbstractWeightedGraph<integer,integer>
{
    view_or_value<std::pair<integer,integer>> adjacentNodes(const integer &u, bool direction) const override
    {
        if(u==1) return std::vector<std::pair<integer,integer >>{};
        if(u%2==0) return std::vector<std::pair<integer,integer>>{{u/2,1}};
        return std::vector<std::pair<integer,integer>>{{3*u+1,1}};
    }

    view_or_value<std::pair<integer,integer>> adjacentNodes(const integer&u) const override
    {
        return adjacentNodes(u,true);
    }

    view_or_value<integer> nodes() const override
    {
        throw std::logic_error("Infinite graph");
    }

    int size() const override
    {
        throw std::logic_error("Infinite graph");
    }
};

int main()
{
    CollatzGraph G;
    integer n;
    std::cin >> n;
    using Set=std::unordered_set<integer>;
    for(auto x:graph::algorithms::reachable_nodes<Set>(G,n))
        std::cout << x << " ";
}
