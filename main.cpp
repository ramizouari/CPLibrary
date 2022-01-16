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
		if (p->v == order_closure<T>{inf}&& p->data.data.has_value() || p->data.data.has_value() && p->v.index()==1 && p->data.data != std::get<T>(p->v) ||
			p->v.index()==1 && !p->data.data.has_value())
			std::cout << std::get<T>(p->v) << ' ';
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
			A.push_back(std::get<T>(p->v));
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
	key_sum_node_t<int,plus_t>* T1 = nullptr, * T2 = nullptr;
	std::cout << "Class Size: " << sizeof(key_sum_node_t<int, plus_t>) << '\n';
	T1 = new key_sum_node_t<int, plus_t>{1,{},nullptr};
	T1->left = new key_sum_node_t<int, plus_t>{ 0,{},T1 };

	T2 = insert(T2, 2);
	T1->h = 2;
	auto T = merge_with_root(T2, T1,(decltype(T1))nullptr);
	b_node<int, int,1000>* b_tree=nullptr;
	constexpr int limit = 1e6;
	std::vector<int> D;
	for (int i = 0; i < limit; i++)
		b_tree = insert(b_tree, i,i);
	to_vector(b_tree, D);
	std::cout << std::endl << "Tree Height: " << height(b_tree) << "\tTree Size: " << size(b_tree);
	for (int i = 0; i < D.size(); i++) if (i != D[i])
		std::cerr << "Error";
	print(b_tree);
	destroy(b_tree);
	_sleep(4000);
}

