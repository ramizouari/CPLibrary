#ifndef __B_TREE_H__
#define __B_TREE_H__
#include "statistic_tree.h"
#include "order.h"
#include <optional>

template<typename T, typename V, int m = 3>
struct b_node;

template<typename T,typename V,int m>
struct b_data
{
	std::optional<V> data;
	b_node<T, V, m>* child;
	b_data() :child(nullptr) {}
	b_data(const V& _data, b_node<T, V, m>* _child = nullptr) :child(nullptr), data(_data) {}
};

template<typename T,typename V,int m>
struct b_node
{
	order_node<order_closure<T>, b_data<T,V,m>>* children;
};

template<typename T, typename V, int m=3>
b_node<T, V, m>* init_b_node()
{
	b_node<T, V, m>* node = new b_node<T, V, m>;
	node->children = nullptr;
	node->children = insert(node->children, inf, b_data<T,V,m>{});
	return node;
}

template<typename T, typename V, int m>
b_node<T, V, m>* rotate_left(b_node<T, V, m>* node, int pos)
{
	auto L = select(node->children, pos-1),R=select(node->children,pos);
	auto s = extract(R, R->v);
	L = insert(L, s->v, s->data);
	delete s;
	return node;
}

template<typename T, typename V, int m>
b_node<T, V, m>* rotate_right(b_node<T, V, m>* node, int pos)
{
	auto L = select(node->children, pos), R = select(node->children, pos+1);
	auto s = extract(L, L->v);
	R = insert(R, s->v, s->data);
	delete s;
	return node;
}

template<typename T, typename V, int m>
b_node<T, V, m>* split(b_node<T, V, m>* node, int pos)
{
	auto A = select(node->children, pos);
	auto Q = extract(A, A->v);
	node->children = Q.second;
	A = Q.first;
	auto [L, R] = split(A, select(A, tree_size(A) / 2)->v);
	return node;
}

template<typename T, typename V, int m>
b_node<T, V, m>* insert_no_root(b_node<T, V, m>* node, const T& v, const V& data)
{

}

template<typename T,typename V,int m>
b_node<T, V, m>* insert(b_node<T, V, m>* node,const T&v,const V& data)
{
	if (!node)
		node=init_b_node<T, V, m>();
	node->children = insert(node->children, v, b_data<T, V, m>{ data,nullptr });
	return node;
}


#endif