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
#include "b_tree.h"

template<typename T,typename V,int m>
void print(b_node<T, V, m>* node)
{
	if (!node)
		return;
	auto p = begin(node->children);
	do
	{
		print(p->data.ptr);
		if (p->v.index()==1)
			std::cout << std::get<T>(p->v) << "->" << p->data.data.value() << '\t';
		std::cout.flush();
	} while (p = next(p));
}

template<typename R,typename T, typename V, int m>
void to_vector(b_node<T, V, m>* node, std::vector<R> &A)
{
	if (!node)
		return;
	auto p = begin(node->children);
	do
	{
		to_vector(p->data.ptr,A);
		if (p->v != order_closure<T>{inf})
			A.push_back(p->data.data.value());
	} while (p = next(p));
}

template<typename T, typename V, int m>
int height(b_node<T, V, m>* node)
{
	if (!node)
		return 0;
	auto p = begin(node->children);
	int h = 0;
	do
	{
		h = std::max(h, height(p->data.ptr) + 1);
	} while (p = next(p));
	return h;
}

template<typename T, typename V, int m>
int size(b_node<T, V, m>* node)
{
	if (!node)
		return 0;
	auto p = begin(node->children);
	int s = 0;
	do
	{
		s+=size(p->data.ptr);
		if (p->v != order_closure<T>{inf})
			s++;
	} while (p = next(p));
	return s;
}

template<typename T,typename V,int m>
void print(order_node<order_closure<T>, b_data<T, V, m>>* node)
{
	if (!node)
		return;
	print(node->left);
	print(node->data.ptr);
	if (node->v != order_closure<T>{inf} && node->data.data.has_value() || node->data.data!=std::get<T>(node->v))
		std::cout << std::get<T>(node->v) << ' ';
	print(node->right);
}

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
	b_node<int, int,10>* b_tree=nullptr;
	constexpr int limit = 100;
	std::vector<int> D;
	for (int i = 0; i < limit; i++)
		b_tree = insert(b_tree, i,i);
	to_vector(b_tree, D);
	std::cout << std::endl << "Tree Height: " << height(b_tree) << "\tTree Size: " << size(b_tree) << '\n';
	for (int i = 0; i < D.size(); i++) if (i != D[i])
		std::cout << "Error: " << i << ' ' << D[i] << '\n';
	for(int i=0;i<limit;i++)
		b_tree = erase(b_tree, i);
	print(b_tree);
	destroy(b_tree);
	_sleep(12000);
}

