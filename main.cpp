#include <iostream>
#include <vector>
#include "nt/modular_arithmetic.h"
#include "nt/number_theory.h"

int L=1e6;

real eps=1e-8;

template<typename T>
auto sorted(std::vector<T> A)
{
    std::sort(A.begin(),A.end());
    return A;
}

int main()
{
   factoriser F(L);
   int n;
   int D=0;
   for(int i=2;i<=L;i++)
       D=std::max<int>(D,F.divisors_count(i));
   std::cout << "Divisors: " << D << '\n';
   std::cin >> n;

   std::vector<real> E(n+1);
   E[1]=0;
   for(auto r:sorted(F.divisors_list(n))) if(r!=1)
   {
       for(auto d:sorted(F.divisors_list(r))) if(d!=r)
            E[r]+=F.totient(r/d)*E[r/d];
       E[r]/=r;
       E[r]=(E[r]+1)/(1.-static_cast<real>(F.totient(r))/r);
   }
   std::cout << "Expected: " << E[n] << '\n';

}