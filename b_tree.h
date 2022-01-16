#ifndef __B_TREE_H__
#define __B_TREE_H__
#include "statistic_tree.h"
#include "order.h"
#include <optional>

/*
* B-tree of an order m
* It is an (m/2,m)-tree which is used mainly on DBMS to minimise disk access
* It is build on top of an Ordered Statistic Tree to support big values of m
* @Requirements:
* - (C,<=) is a totally ordered set, where C is the order closure of T
* - m is even
* @Notes
* Every B-tree contains an element inf, this element does not have a correspending value.
*/
template<typename T, typename V, int m = 4>
struct b_node;


/*
* b_data represents the value mapped by a key k
* It contains two elements:
* 1. data is the value mapped by k.
* 2. ptr contains a B-subtree of the same order with keys less than k
* @Notes
* - The keys in ptr are not less than the previous key in the current order statistic tree
* - Every B-tree contains an element inf, this element does not have a correspending value. In this case the ptr represents
* the B-subtree having all keys not less than the biggest key in the current order statistic tree
*/
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

/*
* Initialise an empty B-Tree
* @Notes
* An empty B-tree has a unique key inf with no data and no children
*/
template<typename T, typename V, int m=4>
b_node<T, V, m>* init_b_node()
{
	b_node<T, V, m>* node = new b_node<T, V, m>;
	node->children = nullptr;
	node->children = insert(node->children, inf, b_data<T,V,m>{});
	return node;
}


/*
* Rotate Right
* @Requirements
* - The right B-subtree must have strictly more than m/2 elemeents
* - The left B-subtree must have strictly less than m elements
*/
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
	q->v = L->v;
	q->data.data=L->data.data;
	L->v = s->v;
	L->data.data = s->data.data;
	s->data.data.reset();
	L->data.ptr->children = insert(L->data.ptr->children, inf, s->data);
	delete s;
	return node;
}


/*
* Rotate Left
* @Requirements
* - The left B-subtree must have strictly more than m/2 elemeents
* - The right B-subtree must have strictly less than m elements
*/
template<typename T, typename V, int m>
b_node<T, V, m>* rotate_right(b_node<T, V, m>* node, int pos)
{
	auto L = select(node->children, pos), R = select(node->children, pos + 1);
	auto p = L->data.ptr->children;
	while (p->right)
		p = p->right;
	auto [s, tmp] = extract(p, p->v);
	L->data.ptr->children = tmp;
	s->data.data=L->data.data;
	s->v = L->v;
	while (tmp->right)
		tmp = tmp->right;
	L->v = tmp->v;
	L->data.data = tmp->data.data;
	tmp->v = inf;
	tmp->data.data.reset();
	R->data.ptr->children = insert(R->data.ptr->children, s->v, s->data);
	delete s;
	return node;
}

/*
* Merge two B-subtrees into one
* @Requirements
* Each B-subtree must have m/2 children
*/
template<typename T, typename V, int m>
b_node<T, V, m>* merge(b_node<T, V, m>* node, int pos)
{
	auto L = select(node->children, pos), R = select(node->children, pos + 1);
	auto p = L->data.ptr->children;
	while (p->right)
		p = p->right;
	p->v = L->v;
	p->data.data = L->data.data;
	R->data.ptr->children = merge(L->data.ptr->children, R->data.ptr->children);
	erase(L, L->v);
	return node;
}


/*
* Split a B-Tree child into 2 children
* @Requirements
* - child is a member of the Ordered Statistic Tree node->children
* - child has m elements
*/
template<typename T, typename V, int m>
b_node<T, V, m>* split(b_node<T, V, m>* node, order_node<order_closure<T>,b_data<T,V,m>> *child)
{
	auto med = select(child->data.ptr->children,m/2);
	auto [Q,tmp] = extract(med, med->v);
	child->data.ptr->children = tmp;
	med = Q;
	auto [L, R] = split(child->data.ptr->children, med->v);
	b_data<T, V, m> b_info_1, b_info_2,placeholder;
	b_info_1 = child->data;
	b_info_1.ptr->children = R;
	b_info_2.data = med->data.data;
	b_info_2.ptr = new b_node<T,V,m>;
	b_info_2.ptr->children = L;
	placeholder.ptr = med->data.ptr;
	b_info_2.ptr->children = insert(b_info_2.ptr->children, inf, placeholder);
	node->children = insert_or_assign(node->children, child->v, b_info_1);
	node->children = insert_or_assign(node->children, med->v, b_info_2);
	return node;
}

/*
* Inserts (v,data) into the B-tree.
*/
template<typename T, typename V, int m>
void insert_no_root(b_node<T, V, m>* node, const typename std::common_type<T>::type & v,
	const typename std::common_type<V>::type& data,bool or_assign=false)
{
	auto L = lower_bound(node->children, v);
	if (or_assign && L->v == order_closure<T>{v})
	{
		L->data.data = data;
		return;
	}
	else if (!L->data.ptr)
	{
		b_data<T, V, m> info;
		info.data = data;
		info.ptr = nullptr;
		node->children = insert(node->children, v, info);
		return;
	}
	else if (size(L->data.ptr->children) == m)
		node = split(node, L);
	insert_no_root(L->data.ptr, v, data,or_assign);
}

template<typename T,typename V,int m>
b_node<T, V, m>* insert(b_node<T, V, m>* node, const typename std::common_type<T>::type& v,
	const typename std::common_type<V>::type& data,bool or_assign = false)
{
	if (!node)
		node=init_b_node<T, V, m>();
	if (size(node->children) == m)
	{
		auto root = init_b_node<T, V, m>();
		root->children->data.ptr = node;
		node = split(root, root->children);
	}
	insert_no_root(node, v, data,or_assign);
	return node;
}

template<typename T, typename V, int m>
b_node<T, V, m>* insert_or_assign(b_node<T, V, m>* node, const typename std::common_type<T>::type& v,
	const typename std::common_type<V>::type& data)
{
	return insert(node, v, data, true);
}

template<typename T, int m>
b_node<T, std::monostate, m>* insert(b_node<T, std::monostate, m>* node, const typename std::common_type<T>::type& v)
{
	return insert(node, v, {});
}

template<typename T, int m>
b_node<T, std::monostate, m>* insert_or_assign(b_node<T, std::monostate, m>* node, const typename std::common_type<T>::type& v)
{
	return insert(node, v, {},true);
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