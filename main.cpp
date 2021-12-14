//
// Created by ramizouari on 01/12/2021.
//
#include "linear_algebra.h"
#include <iostream>
#include "abstract_algebra.h"
#include "fft.h"
#include "ring_extension.h"
#include <chrono>
#include "optimisation.h"


int main()
{
    factoriser F(3e6);
    fast_fourier<>::set_factoriser(F);
    std::vector<IC> u1(pow(5,9),1),u2(1<<21,1);
    auto t1=std::chrono::high_resolution_clock::now();
    fast_fourier<> FFT(u1.size());
    IC r=0;
    auto w=FFT(u1);
    auto t2=std::chrono::high_resolution_clock::now();
    std::cout << "FFT 1: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms\n";
    fast_fourier_base_2<> FFT2(u2.size());
    auto v=FFT2(u2);
    auto t3=std::chrono::high_resolution_clock::now();
    for(auto s:w)
        r+=s;
    for(auto s:v)
        r+=s;
    std::cout << r << '\n';
    std::cout << "FFT 2: " << std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2).count() << "ms";
}