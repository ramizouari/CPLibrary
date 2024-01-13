#include <iostream>
#include <numeric>
#include <map>
#include <fstream>
#include <chrono>
#include "data_structures/fixed/sparse_array.h"
#include "algebra/binary_operation.h"
#include "nt/number_theory.h"
#include "nt/modular/dynamic.h"
#include "topology/optimisation.h"
#include "polynomial/formal_series.h"
#include "signals/ntt.h"
#include "combinatorics/ogf.h"
#include "combinatorics/egf.h"

constexpr cp::integer M=998'244'353,L=1e6;
using IK=cp::cyclic<M>;
int main()
{
    auto F=std::make_shared<cp::light_factoriser>(L);
    cp::default_factoriser_t::default_factoriser = F;

    int n,k;
    std::cin >> n >> k;
    std::vector<IK> X(n);
    X[1]=1;

    auto t1=std::chrono::high_resolution_clock::now();
    auto factorial=std::make_shared<cp::precomputed_factorial<IK>>(L<<1);
    auto Z=cp::egf<M>(X,factorial,false,factorial);
    auto S1=cp::exp_series<M>(1,k+1,factorial,factorial);
    S1.data()[0]=0;
    S1.data().resize(n);
    auto S2=cp::set(S1);
    auto t2=std::chrono::high_resolution_clock::now();
   for(int i=0;i<n;i++)
       std::cout << (cp::integer)S2.weight(i) << ' ';
   std::cerr << '\n' << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << std::endl;
}