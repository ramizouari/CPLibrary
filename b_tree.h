#ifndef __B_TREE_H__
#define __B_TREE_H__
#include "statistic_tree.h"
#include "order.h"
#include <optional>

/*
* B-tree of an order m
* This tree is used mainly on DBMS to optimise disk access*
* It is build on top of an Ordered Statistic Tree to support big values of m
* @Requirements:
* - (T,<=) is a totally ordered set
* - m is even
*/
template<typename T, typename V, int m = 4>
struct b_node;

template<typename T,typename V,int m>
struct b_data
{
	std::optional<V> data;
	b_node<T, V, m>* ptr;
	b_data() :ptr(nullptr) {}
	b_data(const V& _data, b_node<T, V, m>* _child = nullptr) :ptr(nullptr), data(_data) {}
};

template<typename T,typename V,int m>
struct b_node
{
	order_node<order_closure<T>, b_data<T,V,m>>* children;
};

template<typename T, typename V, int m=4>
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
	auto p = R->data.ptr->children;
	while (p->left)
		p = p->left;
	auto [s,tmp] = extract(p, p->v);

	R->data.ptr->children=tmp;
	auto q=lower_bound(L->data.ptr->children, inf);
	q->v = s->v;
	q->data.data=s->data.data;
	s->data.data.reset();
	L->data.ptr->children = insert(L->data.ptr->children, inf, s->data);
	delete s;
	return node;
}

template<typename T, typename V, int m>
b_node<T, V, m>* rotate_right(b_node<T, V, m>* node, int pos)
{
	auto L = select(node->children, pos), R = select(node->children, pos+1);
	while (L->right)
		L = L->right;
	auto [s,tmp] = extract(L, L->v);
	L = tmp;
	R = insert(R, s->v, s->data);
	delete s;
	return node;
}

template<typename T, typename V, int m>
b_node<T, V, m>* split(b_node<T, V, m>* node, int pos)
{
	auto A = select(node->children, pos);
	auto med = select(A->data.ptr->children,m/2);
	auto [Q,tmp] = extract(med, med->v);
	A->data.ptr->children = tmp;
	med = Q;
	auto [L, R] = split(A->data.ptr->children, med->v);
	b_data<T, V, m> b_info_1, b_info_2,placeholder;
	b_info_1 = A->data;
	b_info_1.ptr->children = R;
	b_info_2.data = med->data.data;
	b_info_2.ptr = new b_node<T,V,m>;
	b_info_2.ptr->children = L;
	placeholder.ptr = med->data.ptr;
	b_info_2.ptr->children = insert(b_info_2.ptr->children, inf, placeholder);
	node->children = insert_or_assign(node->children, A->v, b_info_1);
	node->children = insert_or_assign(node->children, med->v, b_info_2);
	return node;
}

template<typename T, typename V, int m>
void insert_no_root(b_node<T, V, m>* node, const typename std::common_type<T>::type & v,
	const typename std::common_type<V>::type& data)
{
	auto L = lower_bound(node->children, v);
	auto pos = order(node->children, L->v);
	if (!L->data.ptr)
	{
		b_data<T, V, m> info;
		info.data = data;
		info.ptr = nullptr;
		node->children = insert(node->children, v, info);
		return;
	}
	else if (size(L->data.ptr->children) == m)
		node = split(node, pos);
	insert_no_root(L->data.ptr, v, data);
}

template<typename T,typename V,int m>
b_node<T, V, m>* insert(b_node<T, V, m>* node, const typename std::common_type<T>::type& v,
	const typename std::common_type<V>::type& data)
{
	if (!node)
		node=init_b_node<T, V, m>();
	if (size(node->children) == m)
	{
		auto root = init_b_node<T, V, m>();
		root->children->data.ptr = node;
		node = split(root, 0);
	}
	insert_no_root(node, v, data);
	return node;
}

template<typename T, int m>
b_node<T, std::monostate, m>* insert(b_node<T, std::monostate, m>* node, const typename std::common_type<T>::type& v)
{
	return insert(node, v, {});
}
template<typename T, typename V, int m>
void destroy(order_node<order_closure<T>, b_data<T,V,m>>* node)
{
	if (!node)
		return;
	destroy(node->left);
	destroy(node->right);
	destroy(node->data.ptr);
	delete node;
}

template<typename T,typename V,int m>
void destroy(b_node<T, V, m>* node)
{
	if (!node)
		return;
	destroy(node->children);
	delete node;
}
#endif