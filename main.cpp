#include "fft.h"
#include "polynomial/polynomial.h"
#include "topology/analysis.h"
#include "functional/zip.h"
#include "data_structures/statistic_tree.h"
#include "algebra/binary_operation.h"
#include "functional/functional.h"
#include "data_structures/data_structures.h"
#include <iostream>
#include <utility>
#include "algebra/order.h"



int main()
{
	ordered_segment_tree<int, plus_t<int>> S;
	for (int i = 0; i < 100; i++)
		S.insert((51 * i + 15) % 100);
	std::cout << S.index_query(0, 21);
}

