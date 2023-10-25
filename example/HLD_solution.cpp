#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <memory>
#include <tuple>
#include <numeric>
#include <functional>
#include <queue>
constexpr class nullopt_t{}nullopt;
template<typename T>
using VE=std::vector<T>;
template<typename T>
class optional
{
    std::unique_ptr<T> value=nullptr;
public:

    optional(const T& x)
    {
        emplace(x);
    }

    optional(nullopt_t)
    {
        value=nullptr;
    }

    optional():optional(nullopt_t{}){}

    optional(T&&x)
    {
        emplace(std::forward<T>(x));
    }
    template<typename ...V>
    void emplace(V &&... x)
    {
        value = std::make_unique<T>(std::forward<V>(x)...);
    }
    T& operator*()
    {
        return *value.get();
    }
    const T& operator*() const
    {
        return *value.get();
    }
    T* operator->()
    {
        return value.get();
    }
    const T* operator->() const
    {
        return value.get();
    }
    bool has_value() const
    {
        return static_cast<bool>(value);
    }
    operator bool() const
    {
        return has_value();
    }
};
template<typename W>
struct WeightedGraph
{
    int n;
    using AdjT=std::pair<int,W>;
    VE<VE<AdjT>> adjacencyList,reverseList;
public:
    explicit WeightedGraph(int _n):n(_n),adjacencyList(n),reverseList(n){}
    void connect(int a,int b, const W & w)
    {
        adjacencyList[a].emplace_back(b,w);
        reverseList[b].emplace_back(a,w);
    }

};

using natural = uint64_t;
using integer = int64_t;
template<typename T>
struct binary_operation
{
    template<typename H0,typename ...H>
    T operator()(const H0&a,const H&... b) const
    {
        if constexpr (sizeof...(b) == 0)
            return a;
        else return reduce(a,this->operator()(b...));
    }
    virtual T reduce(const T& a, const T& b) const = 0;
};
template<typename T>
struct max_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return std::max(a,b);
    }

    inline static T neutral{0};
};
template<typename T>
struct min_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return min(a,b);
    }

    inline static T neutral{};
};
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
template<typename R,typename O>
struct segment_tree
{
    VE<VE<R>> S;
    VE<R> A;
    int n,h;
    segment_tree(const VE<R> &_A):A(_A)
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
    inline static constexpr O F=O();
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

template<typename T,typename O>
struct sparse_array
{
    inline static constexpr O F=O();
    int n,h;
    VE<VE<T>> S;
public:
    sparse_array(const VE<T>&A):n(bit_ceil(A.size())),h(bit_log(n)),S(h+1)
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
template<typename W>
struct WeightedTree : public WeightedGraph<W>
{
    bool reversed=false;
    VE<int> subtree_size;
    using AdjT=typename WeightedGraph<W>::AdjT;
    struct EdgeIdentifier
    {
        int u,v,index;
        bool operator<(const EdgeIdentifier &rhs) const noexcept
        {
            return std::tie(u,v)<std::tie(rhs.u,rhs.v);
        }
        EdgeIdentifier(int _u,int _v,int _index=0):u(_u),v(_v),index(_index){}
    };
    VE<optional<AdjT>> parent;
    int root;
    std::set<EdgeIdentifier> heavy_edges;
    WeightedTree(int n,int _root):WeightedGraph<W>(n),root(_root),subtree_size(n),parent(n)
    {
    }
    explicit WeightedTree(int n):WeightedTree(n,0){}

    void setParent(int u,int v,W w)
    {
        if(reversed) std::swap(u,v);
        this->connect(u,v,w);
        parent[u].emplace(v,w);
    }

    VE<AdjT> &children(int u)
    {
        return reversed?this->adjacencyList[u] : this->reverseList[u] ;
    }

    void buildStatistics()
    {
        updateRoot();
        updateSize(root);
        updateHeavyEdges(root);
        buildLCA();
    }

    void reRoot(int new_root)
    {
        std::queue<int> Q;
        VE<bool> visited(WeightedGraph<W>::n);
        Q.emplace(new_root);
        visited[new_root]=true;
        VE<VE<AdjT>> newAdjacencyList(WeightedGraph<W>::n),newReverseList(WeightedGraph<W>::n);
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
                if(subtree_size[u]>subtree_size[v]) u=parent[u]->first;
                else v=parent[v]->first;
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
        for(int i=0;i<WeightedGraph<W>::n;i++)
            parent[i]=this->adjacencyList[i].empty()?nullopt:optional<AdjT>(this->adjacencyList[i][0]);
        for(int i=0;i<WeightedGraph<W>::n;i++)
            if(!parent[i])
                root=i;
    }

    void updateHeavyEdges(int u)
    {
        for(int i=0;i<children(u).size();i++)
        {
            auto [v,_]=children(u)[i];
            if(subtree_size[v]>=(subtree_size[u]+1)/2)
                heavy_edges.emplace(u,v,i);
            updateHeavyEdges(v);
        }
    }
    void buildLCA()
    {
        VE<HeightData> A;
        euler_tour_endpoints.resize(WeightedGraph<W>::n);
        eulerTour(root,0,A);
        min_t<HeightData>::neutral.first=10000000;
        lca_data=std::make_unique<sparse_array<HeightData,min_t<HeightData>>>(A);
    }

protected:
    using HeightData=std::pair<int,int>;
    using EnpointsData = std::pair<int,int>;
    std::unique_ptr<sparse_array<HeightData,min_t<HeightData>>> lca_data;
    VE<EnpointsData> euler_tour_endpoints;
    void eulerTour(int u,int height,VE<HeightData> &A)
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

template<typename W,typename O, typename RMQ>
struct HeavyLightTree : public WeightedTree<W>
{
private:
    struct RMQIndex
    {
        int rmq_id;
        int index;
        RMQIndex(int _rmq_id,int _index):rmq_id(_rmq_id),index(_index){}
    };
public:
    using EdgeIdentifier=typename WeightedTree<W>::EdgeIdentifier;
    using WeightedTree<W>::WeightedTree;
    using AdjT=typename WeightedTree<W>::AdjT;
    using WeightedTree<W>::root;
    using WeightedTree<W>::parent;
    inline static constexpr O F{};
    void buildStatistics()
    {
        WeightedTree<W>::buildStatistics();
        buildHeavyLightDecomposition();
    }
    W query(int u,int v)
    {
        auto lca=WeightedTree<W>::leastCommonAncestor(u,v);
        return F(query_with_lca(u,lca,false),query_with_lca(v,lca,true));
    }
    void update(int u,W w)
    {
        if(!parent[u])
            throw std::invalid_argument("u must not be the root");
        parent[u]->second=w;
        auto v=parent[u]->first;
        if(WeightedTree<W>::heavy_edges.count(EdgeIdentifier(v,u)))
        {
            auto &R_left=S_left[RMQ_mapper[v].rmq_id];
            auto &R_right=S_right[RMQ_mapper[v].rmq_id];
            R_left.update(component_size[RMQ_mapper[v].rmq_id]-RMQ_mapper[v].index-1,w);
            R_right.update(RMQ_mapper[v].index,w);
        }
    }

    void update(int u,int v, W w)
    {
        if(parent[u] && parent[u]->first==v)
            update(u,w);
        if(parent[v] && parent[v]->first==u)
            update(v,w);
        throw std::invalid_argument("u and v must be adjacent");
    }

private:

    W query_with_lca(int u,int lca, bool invert)
    {
        W R=O::neutral;
        while(RMQ_mapper[u].rmq_id!=RMQ_mapper[lca].rmq_id)
        {
            // x is the upper endpoint of the heavy path containing u
            auto [x,y]=heavy_path_endpoints[RMQ_mapper[u].rmq_id];
            auto b=RMQ_mapper[u].index;
            auto a= RMQ_mapper[x].index;
            auto &S=S_right;
            if(!invert)
                R=F(R,S[RMQ_mapper[u].rmq_id].query(a,b));
            else
                R=F(S[RMQ_mapper[u].rmq_id].query(a,b),R);
            u=x;
            while(u!=lca && !WeightedTree<W>::heavy_edges.count(EdgeIdentifier(parent[u]->first,u)))
            {
                if(!invert) R=F(R,parent[u]->second);
                else R=F(parent[u]->second,R);
                u=parent[u]->first;
            }
        }
        auto b=RMQ_mapper[u].index;
        auto a= RMQ_mapper[lca].index;
        auto &S=S_right;
        if(!invert)
            R=F(R,S[RMQ_mapper[u].rmq_id].query(a,b));
        else
            R=F(S[RMQ_mapper[u].rmq_id].query(a,b),R);
        return R;
    }

    VE<RMQ> S_left, S_right;
    VE<std::pair<int,int>> heavy_path_endpoints;
    VE<int> component_size;
    VE<RMQIndex> RMQ_mapper;
    void buildHeavyLightDecomposition()
    {
        VE<int> stack;
        VE<VE<W>> components;
        stack.push_back(root);
        RMQ_mapper.resize(WeightedTree<W>::n,RMQIndex(-1,-1));
        RMQ_mapper[root]={0,0};
        components.emplace_back();
        heavy_path_endpoints.emplace_back(root,root);
        while(!stack.empty())
        {
            auto u=stack.back();
            stack.pop_back();
            for(auto [v,w]:WeightedTree<W>::children(u))
            {
                if(WeightedTree<W>::heavy_edges.count(EdgeIdentifier(u,v)))
                {
                    auto &C=components[RMQ_mapper[u].rmq_id];
                    auto &[_,y] = heavy_path_endpoints[RMQ_mapper[u].rmq_id];
                    stack.push_back(v);
                    RMQ_mapper[v]={RMQ_mapper[u].rmq_id, RMQ_mapper[u].index+1};
                    y=v;
                    C.push_back(w);
                }
                else
                {
                    heavy_path_endpoints.emplace_back(v,v);
                    stack.push_back(v);
                    RMQ_mapper[v]=RMQIndex(components.size(),0);
                    components.emplace_back();
                }
            }
        }
        for(auto &C:components)
        {
            component_size.push_back(C.size());
            if(C.empty()) C.emplace_back(O::neutral);
            S_right.emplace_back(C);
            reverse(C.begin(),C.end());
            S_left.emplace_back(std::move(C));
        }
    }
};

int main()
{
    using std::cin,std::cout,std::endl;
    int T;
    cin >> T;
    for(int t=1;t<=T;t++)
    {
        int n;
        cin >> n;
        HeavyLightTree<int,max_t<int>,segment_tree<int,max_t<int>>> H(n);
        std::vector<std::pair<int,int>> edges;
        for(int i=0;i<n-1;i++)
        {
            int u,v,w;
            cin >> u >> v >> w;
            u--,v--;
            H.connect(u,v,w);
            edges.emplace_back(u,v);
        }
        H.reRoot(0);
        H.buildStatistics();
        std::string query;
        cin >> query;
        while(query!="DONE")
        {
            if(query=="QUERY")
            {
                int u,v;
                cin >> u >> v;
                u--,v--;
                cout << H.query(u,v) << endl;
            }
            else if(query=="CHANGE")
            {
                int t,w;
                cin >> t >> w;
                t--;
                auto [u,v]=edges[t];
                H.update(u,v,w);
            }
            cin >> query;
        }
    }
}