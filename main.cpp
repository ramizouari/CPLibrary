#include "fft.h"
#include "polynomial.h"
#include "analysis.h"
#include "zip.h"
#include <iostream>
#include <utility>

int main()
{
	s_matrix<IC, 18, 18> A;
	s_vector<int, 18> w;
	s_vector<IC,18> v;
	s_vector<int, 3> p;
	auto&& [a, b, c] = p;
	for (int i = 0; i < 18; i++) for (int j = 0; j < 18; j++)
		A[i][j] = i + j;
	for (int i = 0; i < 18; i++)
	{
		w[i] = 18 - i;
		v[i] = i * i;
	}
	apply_pointwise(std::plus<IC>(), A, A, A);
	std::sort(w.begin(), w.end());
	for (auto [a,b,c,d] : zip(w, v,w,v))
		std::cout << a << ' ' << b << ' ' << c <<  ' ' << d << '\n';

}