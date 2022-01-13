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
#include "ml.h"
#include <fstream>
#include <sstream>
#include <set>

int main()
{
	std::ifstream dataset("iris.data");
	std::vector<std::vector<real>> data;
	std::vector<real> flower_class;
	std::map<std::string, int> encoder;
	int C = 0;
	for (std::string S; std::getline(dataset, S) && S!="";)
	{
		std::stringstream stream(S);
		std::string tmp;
		real a,b,c,d;
		char* send;
		std::getline(stream, tmp, ',');
		a = std::strtold(tmp.c_str(),&send);
		std::getline(stream, tmp, ',');
		b = std::strtold(tmp.c_str(), &send);
		std::getline(stream, tmp, ',');
		c = std::strtold(tmp.c_str(), &send);
		std::getline(stream, tmp, ',');
		d = std::strtold(tmp.c_str(), &send);
		std::string R;
		std::getline(stream, R, ',');
		if (encoder.count(R))
			flower_class.push_back(encoder[R]);
		else
			flower_class.push_back(encoder[R]=C++);
		data.push_back({ a,b,c,d });
	}
	d_matrix<real> X(data);
	d_vector<real> y(flower_class);
	multilogistic_regression M;
	k_nearest_neighbour_classifier<L1_norm<d_vector<real>>> KN;
	auto y_pred= M.fit(X, y).predict(X);
	KN.fit(X, y);
	L2_inner_product<real, d_vector<real>> L2;
	std::cout << '\n' << M.score(X, y) << '\t' << M.error(X,y) << '\n';
	KN.k = 15;
	std::cout << '\n' << KN.score(X, y) << '\n';

	key_sum_node_t<int,plus_t>* T1 = nullptr, * T2 = nullptr;
	std::cout << "Class Size: " << sizeof(key_sum_node_t<int, plus_t>);
	for (int i = 0; i < 1000000; i++)
		T1 = insert(T1, i);
	auto [S1, S2] = split(T1, 500000);

}

