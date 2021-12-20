#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>

using integer = std::int32_t;
using couple = std::pair<integer, integer>;

unsigned int bit_log(unsigned int n)
{
    unsigned char a = 0, b = 30, r = 0;
    while (a <= b)
    {
        auto c = (a + b) / 2;
        if (n >> c)
            a = c + 1;
        else
        {
            b = c - 1;
            r = c - 1;
        }
    }
    if ((1 << (r - 1)) == n)
        return r - 1;
    return r;
}

unsigned int bit_floor(unsigned int n)
{
    return 1 << bit_log(n);
}

unsigned int bit_ceil(unsigned int n)
{
    return 1 << (bit_log(n) + 1);
}

class factoriser
{
    int n;
    std::vector<integer> p_list, smallest_d, p_class, p_order;
public:
    factoriser(int _n) :n(_n), smallest_d(n + 1), p_order(n + 1)
    {
        p_list.reserve(n / log(n));
        std::vector<bool> is_prime(n + 1, true);
        int o = 0;
        integer L = std::ceil(std::sqrt(n));
        smallest_d[2] = 2;
        p_order[2] = o++;
        p_list.push_back(2);
        for (int j = 4; j <= n; j += 2)
        {
            smallest_d[j] = 2;
            is_prime[j] = false;
        }
        smallest_d[3] = 3;
        p_list.push_back(3);
        p_order[3] = o++;
        for (int j = 3; j <= n; j += 3)
        {
            smallest_d[j] = 3;
            is_prime[j] = false;
        }
        int s = 2;
        for (integer i = 5; i <= n; i += s, s = 6 - s) if (is_prime[i])
        {
            p_list.push_back(i);
            smallest_d[i] = i;
            p_order[i] = o++;
            if (i <= L) for (integer j = i * i; j <= n; j += i) if (is_prime[j])
            {
                smallest_d[j] = i;
                is_prime[j] = false;
            }
        }
        p_class.reserve(p_list.size());
        for (auto s : p_list)
        {
            integer C = 1;
            int m = s - 1;
            int r = 0, p = 0;
            while (m > 1)
            {
                if (p != smallest_d[m])
                    C = std::max(C, 1 + p_class[p_order[smallest_d[m]]]);
                p = smallest_d[m];
                m /= smallest_d[m];
            }
            if (p <= 3)
                C = 1;
            p_class.push_back(C);
        }
    }

    int count_prime_factors(integer m) const
    {
        int r = 0, p = 0;
        while (m > 1)
        {
            if (p != smallest_d[m])
                r++;
            p = smallest_d[m];
            m /= smallest_d[m];
        }
        return r;
    }

    const std::vector<integer>& prime_factors(integer m) const
    {
        static std::vector<integer> holder;
        auto p = smallest_divisor(m);
        holder = {};
        while (m > 1)
        {
            auto p = smallest_divisor(m);
            if (holder.empty() || holder.back() != p)
                holder.push_back(p);
            m /= p;
        }
        return holder;
    }

    integer smallest_divisor(integer m) const
    {
        if (m <= n)
            return smallest_d[m];
        integer L = std::ceil(std::sqrt(m));
        for (auto p : p_list)
        {
            if (p > L)
                break;
            else if (m % p == 0)
                return p;
        }
        return m;
    }

    bool is_prime(int m) const
    {
        if (m < n)
            return m > 1 && smallest_d[m] == m;
        else
        {
            integer L = std::ceil(std::sqrt(m));
            for (auto p : p_list)
                if (m % p == 0)
                    return false;
                else if (p > L)
                    break;
            return true;
        }
    }

    auto count_primes() const
    {
        return p_list.size();
    }

    auto prime_lower_bound(int n) const
    {
        return *std::lower_bound(p_list.begin(), p_list.end(), n);
    }

    auto prime_upper_bound(int n) const
    {
        return *std::upper_bound(p_list.begin(), p_list.end(), n);
    }

    int prime_order(int p) const
    {
        return p_order[p];
    }

    auto prime_class(int p) const
    {
        return p_class[prime_order(p)];
    }

    const std::vector<integer>& prime_classes() const
    {
        return p_class;
    }

    const std::vector<integer>& prime_list() const
    {
        return p_list;
    }
};

template<typename R, typename O>
struct segment_tree
{
    std::vector<std::vector<R>> S;
    std::vector<R> A;
    int n, h;
    segment_tree(const std::vector<R>& _A) :A(_A)
    {
        n = bit_ceil(A.size());
        A.resize(n, O::neutral);
        int m = n;
        h = 0;
        while (m)
        {
            m /= 2;
            h++;
        }
        S.resize(h);
        for (int i = 0; i < h; i++)
            S[i].resize(1 << i);
        build();
    }

    segment_tree(std::vector<R>&& _A) :A(std::move(_A))
    {
        n = bit_ceil(A.size());
        A.resize(n, O::neutral);
        int m = n;
        h = 0;
        while (m)
        {
            m /= 2;
            h++;
        }
        S.resize(h);
        for (int i = 0; i < h; i++)
            S[i].resize(1 << i);
        build();
    }

    void update(int i, R u)
    {
        A[i] = u;
        S[h - 1][i] = u;
        int m = h - 2;
        i /= 2;
        while (m >= 0)
        {
            S[m][i] = F(S[m + 1][2 * i], S[m + 1][2 * i + 1]);
            m--;
            i /= 2;
        }
    }

    R query(int l, int r)
    {
        return query(std::max(l, 0), std::min(r, n), 0, n, 0);
    }
private:
    inline static constexpr O F = O();
    void build()
    {
        for (int i = 0; i < n; i++)
            S.back()[i] = A[i];
        for (int i = h - 2; i >= 0; i--) for (int k = 0; k < (1 << i); k++)
            S[i][k] = F(S[i + 1][2 * k], S[i + 1][2 * k + 1]);
    }
    R query(int l, int r, int a, int b, int depth)
    {
        if (l >= r)
            return O::neutral;
        if (l == a && r == b)
            return S[depth][l >> (h - 1 - depth)];
        int mid = (a + b) / 2;
        if (mid > r)
            return query(l, r, a, mid, depth + 1);
        else if (mid < l)
            return query(l, r, mid, b, depth + 1);
        else
            return F(query(l, mid, a, mid, depth + 1), query(mid, r, mid, b, depth + 1));
    }
};



template<typename T, typename V, typename S>
struct statistic_node;

template<typename T, typename V, typename S>
struct statistic_node
{
    T v;
    V data;
    statistic_node* left, * right, * parent;
    int h;
    S statistic;
    statistic_node(T _v, V _data, statistic_node* _parent = nullptr) :v(_v), data(_data), left(nullptr), right(nullptr), parent(_parent), h(1), statistic(v, data) {}
    void update()
    {
        h = std::max(left ? left->h : 0, right ? right->h : 0) + 1;
        S::update(this);
    }
};

template<typename T, typename V, typename S>
int height(statistic_node<T, V, S>* node)
{
    return node ? node->h : 0;
}

template<typename T, typename V, typename S>
int balance(statistic_node<T, V, S>* tree)
{
    return height(tree->left) - height(tree->right);
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* upper_bound(statistic_node<T, V, S>* tree, T v)
{
    if (tree == nullptr)
        return nullptr;
    if (tree->v <= v)
        return upper_bound(tree->right, v);
    else
    {
        auto w = upper_bound(tree->left, v);
        if (w == nullptr)
            return tree;
        return w;
    }
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* reverse_upper_bound(statistic_node<T, V, S>* tree, T v)
{
    if (tree == nullptr)
        return nullptr;
    if (tree->v >= v)
        return reverse_upper_bound(tree->left, v);
    else
    {
        auto w = reverse_upper_bound(tree->right, v);
        if (w == nullptr)
            return tree;
        return w;
    }
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* lower_bound(statistic_node<T, V, S>* tree, T v)
{
    if (tree == nullptr)
        return nullptr;
    if (tree->v < v)
        return lower_bound(tree->right, v);
    else
    {
        auto w = lower_bound(tree->left, v);
        if (w == nullptr)
            return tree;
        return w;
    }
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* reverse_lower_bound(statistic_node<T, V, S>* tree, T v)
{
    if (tree == nullptr)
        return nullptr;
    if (tree->v > v)
        return reverse_lower_bound(tree->left, v);
    else
    {
        auto w = reverse_lower_bound(tree->right, v);
        if (w == nullptr)
            return tree;
        return w;
    }
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* rebalance_right(statistic_node<T, V, S>* x)
{
    auto y = x->left, B = y->right;
    y->right = x;
    y->parent = x->parent;
    if (x->parent)
    {
        if (x->parent->left == x)
            x->parent->left = y;
        else x->parent->right = y;
    }
    x->parent = y;
    x->left = B;
    if (B) B->parent = x;
    x->update();
    return y;
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* rebalance_left(statistic_node<T, V, S>* x)
{
    auto y = x->right, B = y->left;
    y->left = x;
    y->parent = x->parent;
    if (x->parent)
    {
        if (x->parent->left == x)
            x->parent->left = y;
        else x->parent->right = y;
    }
    x->parent = y;
    x->right = B;
    if (B) B->parent = x;
    x->update();
    return y;
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* rebalance(statistic_node<T, V, S>* x)
{
    if (!x)
        return nullptr;
    if (balance(x) < -1)
    {
        if (balance(x->right) == 1)
            rebalance_right(x->right);
        x = rebalance_left(x);
    }
    else if (balance(x) > 1)
    {
        if (balance(x->left) == -1)
            rebalance_left(x->left);
        x = rebalance_right(x);
    }
    x->update();
    if (!x->parent)
        return x;
    return rebalance(x->parent);
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* insert(statistic_node<T, V, S>* tree, T v, V data)
{
    if (!tree)
    {
        tree = new statistic_node<T, V, S>(v, data);
        return tree;
    }
    auto p = lower_bound(tree, v);
    if (p == nullptr)
    {
        p = tree;
        while (p->right) p = p->right;
    }
    else if (p->left)
    {
        p = p->left;
        while (p->right)
            p = p->right;
    }
    auto u = new statistic_node<T, V, S>(v, data, p);
    if (v <= p->v)
        p->left = u;
    else p->right = u;
    p->update();
    return rebalance(p);
}

template<typename T, typename V, typename S>
std::pair<statistic_node<T, V, S>*, statistic_node<T, V, S>*> extract(statistic_node<T, V, S>* tree, T v)
{
    auto p = lower_bound(tree, v);
    if (!p)
        return { nullptr,tree };
    if (!p->left)
    {
        auto w = p->parent ? p->parent : p->right;
        if (p->parent)
        {
            if (p->parent->left == p) p->parent->left = p->right;
            else p->parent->right = p->right;
        }
        if (p->right) p->right->parent = p->parent;
        p->right = nullptr;
        p->parent = nullptr;
        return { p,rebalance(w) };
    }
    else if (!p->left->right)
    {
        auto w = p->left;
        if (p->parent)
        {
            if (p->parent->left == p) p->parent->left = p->left;
            else p->parent->right = p->left;
        }
        if (p->right) p->right->parent = w;
        w->right = p->right;
        w->parent = p->parent;
        p->right = nullptr;
        p->left = nullptr;
        p->parent = nullptr;
        return { p,rebalance(w) };
    }
    else
    {
        auto u = p->left;//Position of replacement
        while (u->right)
            u = u->right;
        auto s = u->parent;//Position of path to be updated
        s->right = u->left;
        if (u->left) u->left->parent = s;
        std::swap(u->v, p->v);
        std::swap(u->data, p->data);
        u->left = nullptr;
        u->right = nullptr;
        u->parent = nullptr;
        return { u,rebalance(s) };
    }

}


template<typename T, typename V, typename S>
statistic_node<T, V, S>* erase(statistic_node<T, V, S>* tree, T v)
{
    auto P = extract(tree, v);
    delete P.first;
    return P.second;
}

template<typename T, typename S>
statistic_node<T, std::tuple<>, S>* insert(statistic_node<T, std::tuple<>, S>* tree, T v)
{
    return insert(tree, v, std::make_tuple());
}


template<typename T, typename V, typename S>
statistic_node<T, V, S>* update(statistic_node<T, V, S>* tree, T v, V data)
{
    auto p = lower_bound(tree, v);
    p->data = data;
    return rebalance(p);
}

template<typename T, typename V, typename S>
void destroy(statistic_node<T, V, S>* node)
{
    if (!node)
        return;
    destroy(node->left);
    destroy(node->right);
    delete node;
}

template<typename V, typename O>
struct sum_stats
{
    inline static constexpr O F = O();
    int size;
    V sum;
    sum_stats() {}
    template<typename T>
    sum_stats(T v, V data) :size(1), sum(data) {}
    template<typename T>
    static void update(statistic_node<T, V, sum_stats>* node);
    template<typename T>
    static int tree_size(statistic_node<T, V, sum_stats>* node);
    template<typename T>
    static V tree_sum(statistic_node<T, V, sum_stats>* node);
};

template<typename V, typename O>
template<typename T>
void sum_stats<V, O>::update(statistic_node<T, V, sum_stats<V, O>>* node) {
    node->statistic.size = (node->left ? node->left->statistic.size : 0) + 1 + (node->right ? node->right->statistic.size : 0);
    node->statistic.sum = F(tree_sum(node->left), F(node->data, tree_sum(node->right)));
}

template<typename V, typename O>
template<typename T>
int sum_stats<V, O>::tree_size(statistic_node<T, V, sum_stats<V, O>>* node)
{
    return node ? node->statistic.size : 0;
}

template<typename V, typename O>
template<typename T>
V sum_stats<V, O>::tree_sum(statistic_node<T, V, sum_stats>* node)
{
    return node ? node->statistic.sum : O::neutral;
}

template<typename T, typename V, typename O>
int order(statistic_node<T, V, sum_stats<V, O>>* tree, T v)
{
    using statistic_type = sum_stats<V, O>;
    if (!tree)
        return 0;
    if (v < tree->v)
        return order(tree->left, v);
    else if (tree->v == v)
    {
        if (tree->right && tree->right->v == v)
            return statistic_type::tree_size(tree->left) + 1 + order(tree->right, v);
        else return statistic_type::tree_size(tree->left);
    }
    else return statistic_type::tree_size(tree->left) + 1 + order(tree->right, v);
}

template<typename T, typename V, typename O>
T select(statistic_node<T, V, sum_stats<V, O>>* tree, int o)
{
    using statistic_type = sum_stats<V, O>;
    int s = statistic_type::tree_size(tree->left);
    if (s == o)
        return tree->v;
    else if (s < o)
        return select(tree->right, o - s - 1);
    else return select(tree->left, o);
}

template<typename T, typename V, typename O>
V prefix_sum(statistic_node<T, V, sum_stats<V, O>>* tree, T U)
{
    using S = sum_stats<V, O>;
    if (!tree)
        return O::neutral;
    else if (tree->v >= U)
        return prefix_sum(tree->left, U);
    else return S::F(S::tree_sum(tree->left), S::F(tree->data, prefix_sum(tree->right, U)));
}

template<typename T, typename V, typename O>
V suffix_sum(statistic_node<T, V, sum_stats<V, O>>* tree, T L)
{
    using S = sum_stats<V, O>;
    if (!tree)
        return O::neutral;
    else if (tree->v < L)
        return suffix_sum(tree->right, L);
    else return S::F(suffix_sum(tree->left, L), S::F(tree->data, S::tree_sum(tree->right)));
}

template<typename T, typename V, typename O>
V sum(statistic_node<T, V, sum_stats<V, O>>* tree, T L, T R)
{
    using S = sum_stats<V, O>;
    if (!tree)
        return O::neutral;
    if (tree->v < L)
        return sum(tree->right, L, R);
    else if (tree->v >= R)
        return sum(tree->left, L, R);
    else return S::F(suffix_sum(tree->left, L), S::F(tree->data, prefix_sum(tree->right, R)));
}

template<typename T>
struct binary_operation
{
    virtual T operator()(const T& a, const T& b) const = 0;
    template<typename H0, typename ...H>
    T operator()(const H0& a, const H&... b)
    {
        return this->operator()(a, this->operator()(b...));
    }
};

template<typename T>
struct plus_t :public binary_operation<T>
{
    T operator()(const T& a, const T& b) const
    {
        return a + b;
    }

    T inv(const T& a) const
    {
        return -a;
    }
    inline static T neutral = T();
};

class forest
{
    std::array < statistic_node<int, int, sum_stats<int, plus_t<int>>>*, 13> S, last;
public:
    forest()
    {
        for (int i = 0; i < 13; i++)
        {
            S[i] = nullptr;
            last[i] = nullptr;
        }
    }
    ~forest()
    {
        for (int i = 0; i < 13; i++)
            destroy(S[i]);
    }

    void insert(int n, int c)
    {
        S[c] = ::insert(last[c], n, 1);
        if (last[c] == nullptr)
            last[c] = S[c];
        else if (last[c]->right)
            last[c] = last[c]->right;
    }

    int query(int a, int b, int c)
    {
        return sum(S[c],a,b);
    }
    std::array<int, 13> query(int a, int b)
    {
        std::array<int, 13> Q;
        for (int i = 0; i < 13; i++)
            Q[i] = query(a, b, i);
        return Q;
    }
};

template<typename R, int n>
class s_vector
{
    std::array<R, n> u;
public:
    using base_field = R;
    using base_ring = R;
    inline static constexpr int dim()
    {
        return n;
    }

    s_vector()
    {
        for (int i = 0; i < n; i++)
            u[i] = 0;
    }

    s_vector(std::array<R, n>_u) :u(std::move(_u)) {}

    auto& operator[](int k)
    {
        return u[k];
    }

    const auto& operator[](int k) const
    {
        return u[k];
    }

    auto& operator+=(const s_vector& o)
    {
        auto r = std::min(dim(), o.dim());
        for (int i = 0; i < r; i++)
            u[i] += o.u[i];
        return *this;
    }

    auto& operator-=(const s_vector& o)
    {
        auto r = std::min(dim(), o.dim());
        for (int i = 0; i < r; i++)
            u[i] -= o.u[i];
        return *this;
    }

    auto& operator*=(R k)
    {
        for (auto& s : u)
            s *= k;
        return *this;
    }

    auto operator+(const s_vector& o) const
    {
        auto v = *this;
        return v += o;
    }

    auto operator-(const s_vector& o) const
    {
        auto v = *this;
        return v -= o;
    }

    auto operator-() const
    {
        auto v = *this;
        for (auto& s : v.u)
            s = -s;
        return v;
    }

    auto& operator/=(R k)
    {
        for (auto& s : u)
            s /= k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto v = *this;
        return v /= k;
    }

    auto begin()
    {
        return u.begin();
    }

    auto begin() const
    {
        return u.cbegin();
    }

    auto end()
    {
        return u.end();
    }

    auto end() const
    {
        return u.cend();
    }
};

template<typename R, int n>
auto operator*(const R& k, const s_vector<R, n>& u)
{
    auto v = u;
    return v *= k;
}

using E = s_vector<int, 13>;
constexpr int L = 50000017;
constexpr std::array<int, 11> P_limit = { 827809,274423,494281,827809 ,1165183,2051041,3160753,
5771863,7592029,12366037,50000017 };

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int T;
    std::cin >> T;
    int k = 0;
    std::vector<int> A(T), B(T);
    for (int i = 0; i < T; i++)
        std::cin >> A[i] >> B[i];
    int max_b = 0;
    for (auto b : B)
        max_b = std::max(max_b, b);
    factoriser F(*std::lower_bound(P_limit.begin(),P_limit.end(),max_b));
    auto d = F.count_primes();
    const std::vector<integer>& C(F.prime_classes()), & P(F.prime_list());
    forest S;
    for (int i = 0; i < d; i++)
        S.insert(P[i], C[i] - 1);
    while (T--)
    {
        int a=A[k], b=B[k];
        k++;
        std::cout << "Case " << k << ":\n";
        auto p = F.prime_order(F.prime_lower_bound(a)), q = F.prime_order(F.prime_upper_bound(b));
        auto R = S.query(a, b+1);
        for (int i = 0; i < 13; i++)
            std::cout << "There are " << R[i] << " primes in class " << i + 1 << ".\n";
        std::cout << '\n';
    }
}