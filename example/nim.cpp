#include <iostream>
#include "linear_algebra/matrix.h"
#include "linear_algebra/special_polynomials.h"
#include "nt/modular_arithmetic.h"
#include "polynomial/fast_polynomial.h"
#include "signals/ntt.h"
#include "nt/primality.h"
#include <random>
#include <chrono>

constexpr int m=1;
constexpr cp::integer M=998244353;
using IK=cp::cyclic<M>;

struct invert_with_cache : public cp::invertible_operation<IK>
{
    mutable std::unordered_map<IK,IK> cache;
    invert_with_cache()
    {
    }
    IK inv(const IK& x) const override
    {
        if(cache.count(x))
            return cache[x];
        return cache[x]=x.inv();
    }
};

int main(int argc, char**argv)
{
    std::ios::sync_with_stdio(false);
    int n,q;
    std::cin >> n >> q;
    cp::linalg::d_matrix<IK> A(0,cp::linalg::m_shape{n,n});
    std::uniform_int_distribution<cp::integer> d(0,M-1);
    std::mt19937_64 rng(std::random_device{}());

    if(argc==1) for(int i=0;i<n;i++) for(int j=0;j<n;j++)
        std::cin >> static_cast<cp::integer&>(A[i][j]);
    else for(int i=0;i<n;i++) for(int j=0;j<n;j++)
        A[i][j]=d(rng);
    auto F=std::make_shared<cp::factoriser>(1e5);
    cp::default_factoriser_t::default_factoriser=F;
    std::vector<cp::linalg::d_vector<IK>> V(m);
    for(int i=0;i<m;i++)
    {
        V[i]=cp::linalg::d_vector<IK>(cp::linalg::v_shape{n});
        for(int j=0;j<n;j++)
            V[i][j]=d(rng);
    }
    cp::polynomial<IK> P;
    auto t1=std::chrono::high_resolution_clock::now();
    if(n>50)
        P=cp::linalg::minimal_polynomial(A,V[0]);
    else
        P=cp::linalg::faddev_lerrier_characteristic_polynomial(A);
    auto t2=std::chrono::high_resolution_clock::now();
    std::cerr << "dt: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << std::endl;
    std::vector<cp::integer> Q(q);
    for(auto &query:Q)
        if(argc==1)
            std::cin >> query;
        else
            query=d(rng);
    std::vector<IK> R;
    R.reserve(q);
    for(auto query:Q)
        R.push_back(P(query));
    for(auto &r:R) {
        if(n&1)
            r=-r;
        std::cout << static_cast<cp::integer>(r) << '\n';
    }
    auto t3=std::chrono::high_resolution_clock::now();
    std::cerr << "dt: " << std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2).count() << std::endl;
}