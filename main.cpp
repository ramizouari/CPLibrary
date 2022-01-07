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
#include "order.h"

int main()
{
	using I=quadratic_extension<int, 0, -1>;
	order_closure<int> O1 = inf,O2=-inf;
	auto O3 = 4*O1 / O2;
	O3 *= 5;
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
	sum_node<order_closure<int>, double, multiplies_t>* tree = nullptr;
	key_sum_node<int, plus_t>* tree2=nullptr;
	tree2 = insert(tree2, 5);
	tree2 = insert(tree2, 15);
	tree = insert(tree, 2, 6);
	tree = insert(tree, 3, 8);
	tree = insert(tree, 3, 15.5);
	std::cout << sum(tree, -inf, inf) << "\t" << order(tree, 3) << '\n';
	std::cout << key_sum(tree2, 6, 16);
	destroy(tree);
	destroy(tree2);
}

