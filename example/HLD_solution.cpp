#include <bits/stdc++.h>
template<typename T>
using VE=std::vector<T>;
using pii=std::pair<int,int>;
namespace cxx14{ constexpr class nullopt_t{} nullopt; template<typename T> class optional { std::unique_ptr<T> p; public: optional():p(nullptr){} optional(nullopt_t):optional(){} optional(const T& _p):p(std::make_unique<T>(_p)){} optional(T&& _p):p(std::make_unique<T>(std::move(_p))){} optional(const optional<T>& _p):p(_p?p?std::make_unique<T>(*_p):nullptr:nullptr){} optional(optional<T>&& _p) noexcept:p(std::move(_p.p)){} optional<T>& operator=(const optional<T>& _p) { if(_p) p=std::make_unique<T>(*_p); else p=nullptr; return *this; } operator bool() const { return p!=nullptr; } T& operator*() { return *p; } const T& operator*() const { return *p; } T* operator->() { return p.get(); } const T* operator->() const { return p.get(); } template<typename ...Z> void emplace(Z&&... z) { p=std::make_unique<T>(std::forward<Z>(z)...); } };}
using integer = std::int64_t;
template<typename T>
struct binary_operation{ using type=T; template<typename H0,typename ...H> T operator()(const H0&a,const H&... b) const { if constexpr (sizeof...(b) == 0) return a; else return reduce(a,this->operator()(b...)); } virtual T reduce(const T& a, const T& b) const = 0; virtual T neutral_element() const { return T{}; }};
template<typename T>struct plus_t:public binary_operation<T>{ T reduce(const T&a,const T&b) const override { return a+b; } inline static T neutral{};};
template<typename T>struct max_t:public binary_operation<T>{ T e; explicit max_t(T _e):e(_e){} max_t(): max_t(neutral){} T reduce(const T&a,const T&b) const override { return std::max(a,b); } inline static T neutral{}; T neutral_element() const override { return e; }};
template<typename T>struct min_t:public binary_operation<T>{ T e; explicit min_t(T _e):e(_e){} min_t(): min_t(neutral){} T reduce(const T&a,const T&b) const override { return std::min(a,b); } inline static T neutral{}; T neutral_element() const override { return e; }};
inline unsigned int bit_log(unsigned int n){ unsigned char a=0,b=30,r=0; while(a<=b) { auto c=(a+b)/2; if(n>>c) a=c+1; else { b=c-1; r=c-1; } } if(r && (1<<(r-1))==n) return r-1; return r;}inline unsigned int bit_floor(unsigned int n){ return 1<<bit_log(n);}inline unsigned int bit_ceil(unsigned int n){ return 1<<(bit_log(n)+1);}
template<typename O>struct segment_tree{ using R=typename O::type; using type=R; VE<VE<R>> S; VE<R> A; int n,h; segment_tree(const VE<R> &_A):A(_A) { n=bit_ceil(A.size()); A.resize(n,O::neutral); int m=n; h=0; while(m) { m/=2; h++; } S.resize(h); for(int i=0;i<h;i++) S[i].resize(1<<i); build(); } void update(int i,R u) { A[i]=u; S[h-1][i]=u; int m=h-2; i/=2; while(m>=0) { S[m][i]=F(S[m+1][2*i],S[m+1][2*i+1]); m--; i/=2; } } R query(int l,int r) { return query(std::max(l,0),std::min(r,n),0,n,0); }private: inline static O F=O(); void build() { for(int i=0;i<n;i++) S.back()[i]=A[i]; for(int i=h-2;i>=0;i--) for(int k=0;k<(1<<i);k++) S[i][k]=F(S[i+1][2*k],S[i+1][2*k+1]); } R query(int l,int r,int a,int b,int depth) { if(l>=r) return O::neutral; if(l==a && r==b) return S[depth][l>>(h-1-depth)]; int mid=(a+b)/2; if(mid>r) return query(l,r,a,mid,depth+1); else if(mid<l) return query(l,r,mid,b,depth+1); else return F(query(l,mid,a,mid,depth+1),query(mid,r,mid,b,depth+1)); }};
template<typename O>struct sparse_array{ using T=typename O::type; using type = T; inline static O F=O(); int n,h; VE<VE<T>> S;public: sparse_array(const VE<T>&A):n(bit_ceil(A.size())),h(bit_log(n)),S(h+1) { int r=1; for(int i=h;i>=0;i--,r*=2) S[i].resize(n-r+1,O::neutral); for(int i=0;i<A.size();i++) S[h][i]=A[i]; r=1; for(int i=h-1;i>=0;i--,r*=2) for(int j=0;j<=n-2*r;j++) S[i][j]=F(S[i+1][j],S[i+1][j+r]); } T query(int l,int r) const { if(l>=r) return O::neutral; auto d=r-l; auto s=bit_floor(d); auto b=bit_log(s); return F(S[h-b][l],S[h-b][r-s]); }};struct HLDIndex{ int hld_id; int index; HLDIndex(int _hld_id, int _index): hld_id(_hld_id), index(_index){}};
template<typename W>struct HLD_t{ VE<bool> is_heavy; VE<std::pair<int,int>> heavy_path_endpoints; VE<int> component_size; VE<HLDIndex> HLD_mapper; VE<VE<W>> components; VE<VE<int>> component_ids;};
enum class TreeStats{ NONE=0, SIZE=0b00001, HEAVY_EDGES=0b00011, LCA=0b00100, HLD=0b01111, RANGE_QUERIES=0b11111,};bool operator &(TreeStats a,TreeStats b){ return static_cast<int>(a) & static_cast<int>(b);}
template<typename W>struct WG
{
    int n;
    using AT=std::pair<int,W>;
    VE<VE<AT>> adjacencyList,reverseList;
    explicit WG(int _n):n(_n),adjacencyList(n),reverseList(n){}
    void connect(int a,int b, const W & w) { adjacencyList[a].emplace_back(b,w); reverseList[b].emplace_back(a,w); }
};
template<typename W>struct WT : public WG<W>
{
    bool reversed=false;
    VE<int> subtree_size;
    using AT=typename WG<W>::AT;
    VE<cxx14::optional<AT>> parent;
    int root; HLD_t<W> HLD; WT(int n,int _root):WG<W>(n),root(_root),subtree_size(n),parent(n) { HLD.is_heavy.resize(n); }
    explicit WT(int n):WT(n,0){}
    void setParent(int u,int v,W w) { if(reversed) std::swap(u,v); this->connect(u,v,w); parent[u].emplace(v,w); }
    VE<AT> &children(int u) { return reversed?this->adjacencyList[u] : this->reverseList[u] ; }
    void buildStatistics(TreeStats stats=TreeStats::HLD) { updateRoot(); if(stats & TreeStats::SIZE) updateSize(root); if(stats & TreeStats::HEAVY_EDGES) updateHeavyEdges(root); if(stats & TreeStats::LCA) buildLCA(); if(stats & TreeStats::HLD) buildHLD_t(); }
    void reRoot(int new_root) { std::queue<int> Q; VE<bool> visited(WG<W>::n); Q.emplace(new_root); visited[new_root]=true; VE<VE<AT>> newAdjacencyList(WG<W>::n),newReverseList(WG<W>::n); while(!Q.empty()) { auto u=Q.front(); Q.pop(); for(auto [v,w]:this->adjacencyList[u]) if(!visited[v]) { visited[v]=true; Q.emplace(v); newReverseList[u].emplace_back(v,w); newAdjacencyList[v].emplace_back(u,w); } for(auto [v,w]:this->reverseList[u]) if(!visited[v]) { visited[v]=true; Q.emplace(v); newReverseList[u].emplace_back(v,w); newAdjacencyList[v].emplace_back(u,w); } } this->adjacencyList=std::move(newAdjacencyList); this->reverseList=std::move(newReverseList); updateRoot(); }
    int leastCommonAncestor(int u,int v) { if(lca_data) { auto [a,b]=euler_tour_endpoints[u]; auto [c,d]=euler_tour_endpoints[v]; if(a>c) { std::swap(a,c); std::swap(b,d); } if(b<d) return lca_data->query(a,c).second; else return lca_data->query(a,d).second; } else { while(u!=v) { if(subtree_size[u]>subtree_size[v]) u=parent[u]->first; else v=parent[v]->first; } return u; } }
    void updateSize(int u) { subtree_size[u]=1; for(auto [v,_]:children(u)) { updateSize(v); subtree_size[u]+=subtree_size[v]; } }
    void updateRoot() { for(int i=0;i<WG<W>::n;i++) parent[i]=this->adjacencyList[i].empty()?cxx14::optional<AT>{}:this->adjacencyList[i][0]; for(int i=0;i<WG<W>::n;i++) if(!parent[i]) root=i; }
    void updateHeavyEdges(int u) { for(int i=0;i<children(u).size();i++) { auto [v,_]=children(u)[i]; if(subtree_size[v]>=(subtree_size[u]+1)/2) HLD.is_heavy[v]=true; updateHeavyEdges(v); } }
    void buildLCA() { VE<HeightData> A; euler_tour_endpoints.resize(WG<W>::n); eulerTour(root,0,A); min_t<HeightData>::neutral.first=std::numeric_limits<int>::max(); lca_data=std::make_unique<sparse_array<min_t<HeightData>>>(A); }
    void buildHLD_t() { VE<int> stack; stack.push_back(root); HLD.HLD_mapper.resize(WT<W>::n, HLDIndex(-1, -1)); HLD.HLD_mapper[root]={0, 0}; HLD.components.emplace_back(); HLD.component_ids.emplace_back(); HLD.heavy_path_endpoints.emplace_back(root,root); HLD.component_ids.back().push_back(root); while(!stack.empty()) { auto u=stack.back(); stack.pop_back(); for(auto [v,w]:children(u)) { if(HLD.is_heavy[v]) { auto &C=HLD.components[HLD.HLD_mapper[u].hld_id]; auto &[_,y] = HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id]; stack.push_back(v); HLD.HLD_mapper[v]={HLD.HLD_mapper[u].hld_id, HLD.HLD_mapper[u].index + 1}; y=v; C.push_back(w); HLD.component_ids[HLD.HLD_mapper[u].hld_id].push_back(v); } else { HLD.heavy_path_endpoints.emplace_back(v,v); stack.push_back(v); HLD.HLD_mapper[v]=HLDIndex(HLD.components.size(), 0); HLD.components.emplace_back(); HLD.component_ids.emplace_back(); HLD.component_ids.back().push_back(v); } } } }
    int distance(int a,int b) { auto lca=leastCommonAncestor(a,b); return distance_with_lca(a,lca)+distance_with_lca(b,lca); }
    int distance_with_lca(int u,int lca) { int d=0; while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id) { if(WT<W>::HLD.is_heavy[u]) { d+=HLD.HLD_mapper[u].index; u=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id].first; } else { d++; u=parent[u]->first; } } return d+= HLD.HLD_mapper[u].index - HLD.HLD_mapper[lca].index; }
    int kth_node(int a,int b,int k)
    {
        auto lca=leastCommonAncestor(a,b);
        auto d1=distance_with_lca(a,lca);
        auto d2= distance_with_lca(b,lca);
        if(k<d1) return kth_parent(a,k);
        else return kth_parent(b,d2-(k-d1));
    }
    int kth_parent(int a, int k)
    {
        int r=0;
        while(r<k)
        {
            if(HLD.is_heavy[a])
            {
                auto t=std::max(HLD.HLD_mapper[a].index-k,0);
                r+=HLD.HLD_mapper[a].index-t;
                a=HLD.component_ids[HLD.HLD_mapper[a].hld_id][t];
            }
            else
            {
                a=parent[a]->first;
                r++;
            }
        }
        return a;
    }
protected: using HeightData=std::pair<int,int>; using EndpointsData = std::pair<int,int>; std::unique_ptr<sparse_array<min_t<HeightData>>> lca_data; VE<EndpointsData> euler_tour_endpoints; void eulerTour(int u,int height,VE<HeightData> &A) { euler_tour_endpoints[u].first=A.size(); for(auto [v,_]: children(u)) { A.emplace_back(height,u); eulerTour(v,height+1,A); } A.emplace_back(height,u); euler_tour_endpoints[u].second=A.size(); }};
template<typename O, typename RMQ>struct CHLT : public WT<typename O::type>{ using W=typename O::type; using type=W;private:public: using WT<W>::HLD; using WT<W>::WT; using AT=typename WT<W>::AT; using WT<W>::root; using WT<W>::parent; inline static O F{}; void buildStatistics(TreeStats stats=TreeStats::RANGE_QUERIES) { WT<W>::buildStatistics(stats); if(stats & TreeStats::RANGE_QUERIES) buildRangeStatistics(); } W query(int u,int v) { auto lca=WT<W>::leastCommonAncestor(u,v); return F(query_with_lca(u,lca),query_with_lca(v,lca)); } void update(int u,int v, W w) { if(parent[u] && parent[u]->first==v) return update(u,w); if(parent[v] && parent[v]->first==u) return update(v,w); throw std::invalid_argument("u and v must be adjacent"); } void update(int u,W w) { if(!parent[u]) throw std::invalid_argument("u must not be the root"); parent[u]->second=w; auto v=parent[u]->first; if(WT<W>::HLD.is_heavy[u]) { auto &R=S[HLD.HLD_mapper[v].hld_id]; R.update(HLD.HLD_mapper[v].index, w); } }protected: W query_with_lca(int u,int lca) { if(u==lca) return O::neutral; W R=O::neutral; while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id) { auto [x,y]=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id]; auto b=HLD.HLD_mapper[u].index; auto a= HLD.HLD_mapper[x].index; R=F(S[HLD.HLD_mapper[u].hld_id].query(a, b), R); u=x; while(u!=lca && !WT<W>::HLD.is_heavy[u]) { R=F(R,parent[u]->second); u=parent[u]->first; } } auto b=HLD.HLD_mapper[u].index; auto a= HLD.HLD_mapper[lca].index; R=F(R,S[HLD.HLD_mapper[u].hld_id].query(a, b)); return R; } VE<RMQ> S; void buildRangeStatistics() { for(auto &C:HLD.components) { HLD.component_size.push_back(C.size()); if(C.empty()) C.emplace_back(O::neutral); S.emplace_back(C); } }};
int main()
{
    int T;
    using namespace std; std::ios_base::sync_with_stdio(false); cin.tie(nullptr);
    cin >> T;
    for(int t=1;t<=T;t++)
    {
        int n; cin >> n;
        CHLT<plus_t<integer>,segment_tree<plus_t<integer>>> H(n);
        VE<std::pair<int,int>> edges;
        for(int i=0;i<n-1;i++)
        {
            int u,v,w; cin >> u >> v >> w;
            u--,v--; H.connect(u,v,w);
            edges.emplace_back(u,v);
        }
        H.reRoot(0);
        H.buildStatistics();
        string query; cin >> query;
        while(query!="DONE")
        {
            if(query=="DIST") { int u,v; cin >> u >> v; u--,v--; cout << H.query(u,v) << '\n'; }
            else if(query=="KTH") { int u,v,k; cin >> u >> v >> k; u--,v--,k--; cout << H.kth_node(u,v,k)+1 << '\n'; }
            cin >> query;
        }
        cout << '\n';
    }
}