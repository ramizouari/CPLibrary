#include "fft.h"
#include "polynomial.h"
#include "analysis.h"
#include "zip.h"
#include "statistic_tree.h"
#include "binary_operation.h"
#include "functional.h"
#include "data_structures.h"
#include <iostream>
#include <utility>

int main()
{
	using I=quadratic_extension<int, 0, -1>;
	s_matrix<I, 18, 18> A;
	s_matrix<I, 3, 3> B;
	s_vector<I, 18> w;
	s_vector<I,18> v;
	s_vector<I, 3> p;
	auto &[a, b, c] = p;
	a = 5;
	for (int i = 0; i < 18; i++) for (int j = 0; j < 18; j++)
		A[i][j] = i + j;
	for (int i = 0; i < 18; i++)
	{
		w[i] = 18 - i;
		v[i] = i * i;
	}
	apply_pointwise(multiplies_t<I>(), A, A, A);
	statistic_node<int, int, sum_stats<int, multiplies_t<int>>>* tree = nullptr;
	tree = insert(tree, 2, 6);
	tree = insert(tree, 3, 8);
	tree = insert(tree, 3, 15);
	std::cout << sum(tree, 2, 20) << "\t" << order(tree,3);
	destroy(tree);
}

