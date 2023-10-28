#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include "nt/number_theory.h"
constexpr int L=1e6+1;
int main()
{
    factoriser F(L);
    std::vector<integer> A(L);
    integer n;
    std::cin >> n;
    integer R=0;
    for(int i=1;i<L && R < n*n/2;i++)
    {
        A[i]=F.totient(i);
        R+=A[i];
    }
}