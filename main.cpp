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
	ordered_segment_tree<int, plus_t<int>> S;
	for (int i = 0; i < 100; i++)
		S.insert((51 * i + 15) % 100);
	std::cout << S.index_query(0, 21);
}

