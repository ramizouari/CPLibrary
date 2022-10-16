#include <iostream>
#include <vector>
#include <polynomial/ring_extension.h>
#include <nt/modular_arithmetic.h>
#include <polynomial/fft.h>


constexpr integer M=998244353 ;
using IK=d_cyclic;
using R=polynomial<IK>;


void Q_divide(R&p,integer n,integer k,integer d,integer iter)
{

    for(integer i=n;i>=d;i--) for(integer j=k;j<=i;j+=k)
        p[i]+=p[i-j];
}


int main()
{
    factoriser F(1e5);
    fast_ntt<>::set_factoriser(F);
    d_cyclic::m=M;
    R P,w=1;
    integer n,k;
    std::cin >> n >> k;
    P.p.resize(n+1);
    w.p.resize(n+1);
    for(integer i=1,d=i*(2*k+i-1)/2;i<=n-k && d<=n;i++,d=i*(2*k+i-1)/2)
    {
        Q_divide(w,n,k+i-1,d,i);
        for(integer j=0;j+d<=n && j<=w.degree();j++)
            P[j+d]+=w[j];

    }
    integer d = (n-k)*n/2;
    Q_divide(w,n,n-k,d,n-k);
    for(int i=0;i+d<=n && i<=w.degree();i++)
        P[i+d]+=w[i];
    for(int i=1;i<=n;i++)
        std::cout << (integer)P[i] << ' ';
}