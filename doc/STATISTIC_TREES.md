# Statistic Trees Documentation

This document provides comprehensive documentation for the statistic tree data structures implemented in CPLibrary. Statistic trees are augmented AVL trees that maintain additional statistics about subtrees, enabling efficient order statistics, rank queries, and aggregate range operations.

## Table of Contents

- [Overview](#overview)
- [Core Concepts](#core-concepts)
- [Order Statistic Tree](#order-statistic-tree)
- [Key Sum Statistic Tree](#key-sum-statistic-tree)
- [Value Sum Statistic Tree](#value-sum-statistic-tree)
- [B-Tree](#b-tree)
- [Implementation Details](#implementation-details)
- [API Reference](#api-reference)
- [Usage Examples](#usage-examples)
- [Comparison with Other Structures](#comparison-with-other-structures)
- [Advanced Operations](#advanced-operations)

## Overview

### What are Statistic Trees?

**Statistic trees** are self-balancing binary search trees (based on AVL trees) augmented with additional statistical information at each node. This augmentation allows efficient computation of aggregate functions over ranges and support for order-based queries.

**Key Features:**
- **Self-Balancing**: AVL tree rotations maintain O(log n) height
- **Augmented Data**: Each node stores additional statistics about its subtree
- **Generic Statistics**: Support for various statistics (size, sum, etc.)
- **Efficient Queries**: O(log n) for order statistics, rank, and range operations

### Why Use Statistic Trees?

Traditional BSTs support:
- Insert, delete, search in O(log n)

Statistic trees additionally support:
- **Order statistics**: Find kth smallest/largest element in O(log n)
- **Rank queries**: Count elements less than a value in O(log n)
- **Range aggregates**: Compute sum/product over key/value ranges in O(log n)
- **Median queries**: Find median in O(log n)
- **Split/merge**: Split tree at threshold, merge two trees in O(log n)

---

## Core Concepts

### Statistic Node Structure

```cpp
template<typename T, typename V, typename S>
struct statistic_node {
    T v;                    // Key
    V data;                 // Value
    int h;                  // Height (for AVL balancing)
    S statistic;            // Augmented statistic
    statistic_node *left;   // Left child
    statistic_node *right;  // Right child
    statistic_node *parent; // Parent pointer
};
```

**Key Components:**
- **T**: Key type (must be totally ordered)
- **V**: Value type (can be `std::monostate` for key-only trees)
- **S**: Statistic type (defines what to track)

### Statistic Types

CPLibrary provides three main statistic types:

| Statistic | Tracks | Use Cases |
|-----------|--------|-----------|
| **order_stats** | Subtree size | Order statistics, rank, select, median |
| **key_sum_stats** | Aggregate over keys | Range sum/product of keys |
| **sum_stats** | Aggregate over values | Range sum/product of values |

### Statistic Requirements

A statistic type `S` must provide:

```cpp
struct my_statistic {
    int size;  // Required: subtree size

    // Constructor from key and value
    my_statistic(const T& key, const V& value);

    // Static update method
    static void update(statistic_node<T,V,my_statistic>* node);
};
```

The `update` method recomputes statistics after tree modifications.

---

## Order Statistic Tree

### Purpose

Order statistic trees support **order-based queries**:
- Find the kth smallest element
- Determine the rank (position) of an element
- Find median efficiently
- Count elements in a range

### Implementation

**Location**: `include/data_structures/statistic_tree/order_statistics.h`

#### Statistic Structure

```cpp
struct order_stats {
    int size;  // Number of nodes in subtree

    order_stats() {}

    template<typename T, typename V>
    order_stats(const T& v, const V& data) : size(1) {}

    template<typename T, typename V>
    static void update(statistic_node<T,V,order_stats>* node) {
        node->statistic.size =
            (node->left ? node->left->statistic.size : 0) +
            1 +
            (node->right ? node->right->statistic.size : 0);
    }
};
```

**Key Idea**: Each node maintains the size of its subtree, enabling navigation by position.

### Time and Space Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Insert | O(log n) | O(1) |
| Delete | O(log n) | O(1) |
| Find | O(log n) | O(1) |
| Select (kth) | O(log n) | O(1) |
| Rank | O(log n) | O(1) |
| Median | O(log n) | O(1) |
| Update | O(log n) | O(1) |
| Memory per node | - | O(1) |
| Total memory | - | O(n) |

### API Reference

```cpp
namespace cp::data_structures::stats_trees {

// Type alias for order statistic nodes
template<typename T, typename V = std::monostate>
using order_node = statistic_node<T, V, order_stats>;

// Get subtree size
template<typename T, typename V, typename OrderStats>
int size(statistic_node<T, V, OrderStats>* node);

// Select kth smallest element (0-indexed)
template<typename T, typename V, typename OrderStats>
statistic_node<T,V,OrderStats>* select(
    statistic_node<T,V,OrderStats>* tree,
    int k
);

// Get rank (position) of element with key v
template<typename T, typename V, typename OrderStats>
int order(
    statistic_node<T,V,OrderStats>* tree,
    const T& v
);

// Find median element
template<typename T, typename V, typename OrderStats>
statistic_node<T,V,OrderStats>* median(
    statistic_node<T,V,OrderStats>* tree
);

// Rank of smallest element >= v
template<typename T, typename V, typename OrderStats>
int order_inf(
    statistic_node<T,V,OrderStats>* tree,
    const T& v
);

// Rank of largest element <= v
template<typename T, typename V, typename OrderStats>
int order_sup(
    statistic_node<T,V,OrderStats>* tree,
    const T& v
);

}
```

### Algorithm: Select (Find kth Element)

```cpp
template<typename T, typename V, typename OrderStats>
statistic_node<T,V,OrderStats>* select(
    statistic_node<T,V,OrderStats>* tree,
    int k
) {
    int s = size(tree->left);  // Size of left subtree

    if(s == k)
        return tree;  // Current node is kth
    else if(s < k)
        return select(tree->right, k - s - 1);  // Search right
    else
        return select(tree->left, k);  // Search left
}
```

**Intuition**:
- If left subtree has `k` elements, current node is kth
- If left subtree has `< k` elements, recurse right
- Otherwise, recurse left

### Algorithm: Order (Find Rank)

```cpp
template<typename T, typename V, typename OrderStats>
int order(statistic_node<T,V,OrderStats>* tree, const T& v) {
    if(!tree)
        return 0;

    if(v < tree->v)
        return order(tree->left, v);
    else if(tree->v == v) {
        // Handle duplicates
        if(tree->right && tree->right->v == v)
            return size(tree->left) + 1 + order(tree->right, v);
        else
            return size(tree->left);
    }
    else
        return size(tree->left) + 1 + order(tree->right, v);
}
```

### Usage Examples

#### Example 1: Basic Order Statistics

```cpp
#include "data_structures/statistic_tree.h"
#include <iostream>

using namespace cp::data_structures::stats_trees;

int main() {
    order_node<int>* tree = nullptr;

    // Insert elements
    tree = insert(tree, 50, {});
    tree = insert(tree, 30, {});
    tree = insert(tree, 70, {});
    tree = insert(tree, 20, {});
    tree = insert(tree, 40, {});
    tree = insert(tree, 60, {});
    tree = insert(tree, 80, {});

    // Tree contains: {20, 30, 40, 50, 60, 70, 80}

    // Select 3rd smallest (0-indexed)
    auto node = select(tree, 3);
    std::cout << "3rd smallest (0-indexed): " << node->v << std::endl;
    // Output: 50

    // Find rank of 60
    int rank = order(tree, 60);
    std::cout << "Rank of 60: " << rank << std::endl;
    // Output: 4

    // Get median
    auto median_node = median(tree);
    std::cout << "Median: " << median_node->v << std::endl;
    // Output: 50

    // Total size
    std::cout << "Total elements: " << size(tree) << std::endl;
    // Output: 7

    // Cleanup
    destroy(tree);
    return 0;
}
```

#### Example 2: Dynamic Order Queries

```cpp
#include "data_structures/statistic_tree.h"
#include <iostream>

using namespace cp::data_structures::stats_trees;

int main() {
    order_node<int>* tree = nullptr;

    // Dynamic insertions with queries
    std::vector<int> elements = {15, 10, 20, 8, 12, 16, 25};

    for(int i = 0; i < elements.size(); i++) {
        tree = insert(tree, elements[i], {});

        // Query median after each insertion
        auto med = median(tree);
        std::cout << "After inserting " << elements[i]
                  << ", median is " << med->v << std::endl;
    }

    // Output:
    // After inserting 15, median is 15
    // After inserting 10, median is 10
    // After inserting 20, median is 15
    // After inserting 8, median is 10
    // After inserting 12, median is 12
    // After inserting 16, median is 15
    // After inserting 25, median is 15

    destroy(tree);
    return 0;
}
```

#### Example 3: Range Count Queries

```cpp
#include "data_structures/statistic_tree.h"
#include <iostream>

using namespace cp::data_structures::stats_trees;

int main() {
    order_node<int>* tree = nullptr;

    // Insert elements
    int arr[] = {5, 2, 8, 1, 9, 3, 7, 4, 6};
    for(int x : arr) {
        tree = insert(tree, x, {});
    }

    // Count elements in range [3, 7]
    // This is order(8) - order(3)
    int count = order(tree, 8) - order(tree, 3);
    std::cout << "Elements in [3, 8): " << count << std::endl;
    // Output: 5 (elements: 3, 4, 5, 6, 7)

    // Find smallest element >= 5
    auto node = lower_bound(tree, 5);
    std::cout << "Smallest >= 5: " << node->v << std::endl;
    // Output: 5

    // Find largest element < 5
    auto node2 = reverse_upper_bound(tree, 5);
    std::cout << "Largest < 5: " << node2->v << std::endl;
    // Output: 4

    destroy(tree);
    return 0;
}
```

#### Example 4: Key-Value Mapping with Order Statistics

```cpp
#include "data_structures/statistic_tree.h"
#include <iostream>
#include <string>

using namespace cp::data_structures::stats_trees;

int main() {
    // Order statistic tree with string values
    order_node<int, std::string>* tree = nullptr;

    // Insert (key, value) pairs
    tree = insert(tree, 100, std::string("Alice"));
    tree = insert(tree, 85, std::string("Bob"));
    tree = insert(tree, 92, std::string("Charlie"));
    tree = insert(tree, 78, std::string("David"));
    tree = insert(tree, 95, std::string("Eve"));

    // Find student with 3rd highest score
    auto node = select(tree, size(tree) - 3);
    std::cout << "3rd highest score: " << node->v
              << " (" << node->data << ")" << std::endl;
    // Output: 3rd highest score: 92 (Charlie)

    // Find rank of score 92
    int rank = order(tree, 92);
    std::cout << "Rank of 92: " << rank << std::endl;
    // Output: 2

    destroy(tree);
    return 0;
}
```

### Special Features

1. **Duplicate Handling**: Supports duplicate keys with proper rank calculation
2. **Median in O(log n)**: Direct median queries without sorting
3. **Dynamic Queries**: Efficient queries after each insertion/deletion
4. **Range Counting**: Count elements in any range efficiently

---

## Key Sum Statistic Tree

### Purpose

Key sum statistic trees maintain **aggregates over keys** in addition to order statistics. This enables:
- Range sum/product of keys in [L, R)
- Aggregate queries by rank (first k keys)
- Efficient computation of key statistics

### Implementation

**Location**: `include/data_structures/statistic_tree/key_statistics.h`

#### Statistic Structure

```cpp
template<typename T, typename O>
struct key_sum_stats {
    inline static O F = O();  // Binary operation
    int size;
    T key_sum;  // Aggregate of keys in subtree

    template<typename V>
    key_sum_stats(T v, V data) : size(1), key_sum(v) {}

    template<typename V>
    static void update(statistic_node<T, V, key_sum_stats>* node) {
        node->statistic.size =
            (node->left ? node->left->statistic.size : 0) +
            1 +
            (node->right ? node->right->statistic.size : 0);

        node->statistic.key_sum = F(
            tree_key_sum(node->left),
            node->v,
            tree_key_sum(node->right)
        );
    }

    inline static const T& key_neutral = O::neutral;
};
```

**Key Idea**: Maintain aggregate of all keys in subtree using associative operation.

### Time and Space Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Insert | O(log n) | O(1) |
| Delete | O(log n) | O(1) |
| Range Key Sum [L, R) | O(log n) | O(1) |
| Prefix Key Sum [0, R) | O(log n) | O(1) |
| Index Key Sum [first k] | O(log n) | O(1) |
| Memory per node | - | O(1) |

### API Reference

```cpp
namespace cp::data_structures::stats_trees {

// Type aliases
template<typename T, template<typename> typename O, typename V=std::monostate>
using key_sum_node_t = statistic_node<T, V, key_sum_stats<T, O<T>>>;

template<typename T, typename O, typename V = std::monostate>
using key_sum_node = statistic_node<T, V, key_sum_stats<T, O>>;

// Get sum of all keys in subtree
template<typename T, typename V, typename KeySumStats>
T tree_key_sum(statistic_node<T, V, KeySumStats>* node);

// Sum of keys < U
template<typename T, typename V, typename KeySumStats>
T prefix_key_sum(
    statistic_node<T, V, KeySumStats>* tree,
    const T& U
);

// Sum of keys >= L
template<typename T, typename V, typename KeySumStats>
T suffix_key_sum(
    statistic_node<T, V, KeySumStats>* tree,
    const T& L
);

// Sum of keys in [L, R)
template<typename T, typename V, typename KeySumStats>
T key_sum(
    statistic_node<T, V, KeySumStats>* tree,
    const T& L,
    const T& R
);

// Sum of first n keys (by order)
template<typename T, typename V, typename KeySumStats>
T prefix_index_key_sum(
    statistic_node<T, V, KeySumStats>* tree,
    int n
);

// Sum of keys at indices [a, b)
template<typename T, typename V, typename KeySumStats>
T index_key_sum(
    statistic_node<T, V, KeySumStats>* tree,
    int a,
    int b
);

}
```

### Usage Examples

#### Example 1: Range Sum of Keys

```cpp
#include "data_structures/statistic_tree.h"
#include "algebra/binary_operation.h"
#include <iostream>

using namespace cp;
using namespace cp::data_structures::stats_trees;

int main() {
    // Tree with sum statistic on keys
    key_sum_node_t<int, plus_t>* tree = nullptr;

    // Insert keys
    int keys[] = {10, 5, 15, 3, 7, 12, 20};
    for(int k : keys) {
        tree = insert(tree, k, {});
    }

    // Sum of all keys
    int total = tree_key_sum(tree);
    std::cout << "Sum of all keys: " << total << std::endl;
    // Output: 72

    // Sum of keys in range [5, 15)
    int range_sum = key_sum(tree, 5, 15);
    std::cout << "Sum of keys [5, 15): " << range_sum << std::endl;
    // Output: 5 + 7 + 10 + 12 = 34

    // Sum of 3 smallest keys
    int first_3 = prefix_index_key_sum(tree, 3);
    std::cout << "Sum of 3 smallest: " << first_3 << std::endl;
    // Output: 3 + 5 + 7 = 15

    destroy(tree);
    return 0;
}
```

#### Example 2: Product of Keys

```cpp
#include "data_structures/statistic_tree.h"
#include "algebra/binary_operation.h"
#include <iostream>

using namespace cp;
using namespace cp::data_structures::stats_trees;

int main() {
    // Tree with product statistic on keys
    key_sum_node_t<int, multiplies_t>* tree = nullptr;

    // Insert keys
    tree = insert(tree, 2, {});
    tree = insert(tree, 3, {});
    tree = insert(tree, 5, {});
    tree = insert(tree, 7, {});

    // Product of all keys
    int product = tree_key_sum(tree);
    std::cout << "Product of all keys: " << product << std::endl;
    // Output: 2 * 3 * 5 * 7 = 210

    // Product of keys in [3, 7)
    int range_product = key_sum(tree, 3, 7);
    std::cout << "Product [3, 7): " << range_product << std::endl;
    // Output: 3 * 5 = 15

    destroy(tree);
    return 0;
}
```

---

## Value Sum Statistic Tree

### Purpose

Value sum statistic trees maintain **aggregates over values** stored at nodes. This enables:
- Range sum/product of values for keys in [L, R)
- Aggregate queries by rank
- Efficient value statistics

### Implementation

**Location**: `include/data_structures/statistic_tree/value_stats.h`

#### Statistic Structure

```cpp
template<typename V, typename O>
struct sum_stats {
    inline static O F = O();  // Binary operation
    int size;
    V sum;  // Aggregate of values in subtree

    template<typename T>
    sum_stats(T v, V data) : size(1), sum(data) {}

    template<typename T>
    static void update(statistic_node<T, V, sum_stats>* node) {
        node->statistic.size =
            (node->left ? node->left->statistic.size : 0) +
            1 +
            (node->right ? node->right->statistic.size : 0);

        node->statistic.sum = F(
            tree_sum(node->left),
            node->data,
            tree_sum(node->right)
        );
    }

    inline static const V& neutral = O::neutral;
};
```

### Time and Space Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Insert | O(log n) | O(1) |
| Delete | O(log n) | O(1) |
| Range Sum [L, R) | O(log n) | O(1) |
| Index Sum [first k] | O(log n) | O(1) |
| Update Value | O(log n) | O(1) |
| Memory per node | - | O(1) |

### API Reference

```cpp
namespace cp::data_structures::stats_trees {

// Type aliases
template<typename T, typename V, template<typename> typename O>
using sum_node_t = statistic_node<T, V, sum_stats<V, O<V>>>;

template<typename T, typename V, typename O>
using sum_node = statistic_node<T, V, sum_stats<V, O>>;

// Get sum of all values in subtree
template<typename T, typename V, typename SumStats>
V tree_sum(statistic_node<T, V, SumStats>* node);

// Sum of values for keys < U
template<typename T, typename V, typename SumStats>
V prefix_sum(
    statistic_node<T, V, SumStats>* tree,
    const T& U
);

// Sum of values for keys >= L
template<typename T, typename V, typename SumStats>
V suffix_sum(
    statistic_node<T, V, SumStats>* tree,
    const T& L
);

// Sum of values for keys in [L, R)
template<typename T, typename V, typename SumStats>
V sum(
    statistic_node<T, V, SumStats>* tree,
    const T& L,
    const T& R
);

// Sum of first n values (by key order)
template<typename T, typename V, typename SumStats>
V prefix_index_sum(
    statistic_node<T, V, SumStats>* tree,
    int n
);

// Sum of values at indices [a, b)
template<typename T, typename V, typename SumStats>
V index_sum(
    statistic_node<T, V, SumStats>* tree,
    int a,
    int b
);

}
```

### Usage Examples

#### Example 1: Range Sum of Values

```cpp
#include "data_structures/statistic_tree.h"
#include "algebra/binary_operation.h"
#include <iostream>

using namespace cp;
using namespace cp::data_structures::stats_trees;

int main() {
    // Tree mapping scores to counts
    sum_node_t<int, int, plus_t>* tree = nullptr;

    // Insert (score, count) pairs
    tree = insert(tree, 85, 5);   // 5 students scored 85
    tree = insert(tree, 90, 3);   // 3 students scored 90
    tree = insert(tree, 75, 7);   // 7 students scored 75
    tree = insert(tree, 95, 2);   // 2 students scored 95
    tree = insert(tree, 80, 4);   // 4 students scored 80

    // Total students
    int total = tree_sum(tree);
    std::cout << "Total students: " << total << std::endl;
    // Output: 21

    // Students with score in [80, 95)
    int range_count = sum(tree, 80, 95);
    std::cout << "Students in [80, 95): " << range_count << std::endl;
    // Output: 4 + 5 + 3 = 12

    // Students with top 2 scores
    int top2 = suffix_index_sum(tree, 2);
    std::cout << "Top 2 scores count: " << top2 << std::endl;
    // Output: 2 + 3 = 5

    destroy(tree);
    return 0;
}
```

#### Example 2: Weighted Statistics

```cpp
#include "data_structures/statistic_tree.h"
#include "algebra/binary_operation.h"
#include <iostream>

using namespace cp;
using namespace cp::data_structures::stats_trees;

int main() {
    // Tree mapping items to weights
    sum_node_t<int, double, plus_t>* tree = nullptr;

    // Insert (item_id, weight) pairs
    tree = insert(tree, 10, 5.5);
    tree = insert(tree, 20, 3.2);
    tree = insert(tree, 30, 7.8);
    tree = insert(tree, 40, 2.1);
    tree = insert(tree, 50, 6.4);

    // Total weight
    double total_weight = tree_sum(tree);
    std::cout << "Total weight: " << total_weight << std::endl;
    // Output: 25.0

    // Weight of items [20, 40)
    double range_weight = sum(tree, 20, 40);
    std::cout << "Weight [20, 40): " << range_weight << std::endl;
    // Output: 3.2 + 7.8 = 11.0

    destroy(tree);
    return 0;
}
```

---

## B-Tree

### Purpose

B-trees are multi-way search trees optimized for systems with slow disk access. This implementation uses order statistic trees internally to handle large branching factors efficiently.

### Implementation

**Location**: `include/data_structures/b_tree.h`

#### Structure

```cpp
template<typename T, typename V, int m = 4>
struct b_node {
    // Children stored in order statistic tree
    stats_trees::order_node<order_closure<T>, b_data<T,V,m>>* children;
};

template<typename T, typename V, int m>
struct b_data {
    std::optional<V> data;  // Value (optional for sentinel)
    b_node<T, V, m>* ptr;   // Child pointer
};
```

**Key Idea**: Use order statistic trees to efficiently manage the (potentially large) number of children in each B-tree node.

### Time Complexity

| Operation | Time |
|-----------|------|
| Insert | O(log n) |
| Delete | O(log n) |
| Find | O(log n) |
| Range Query | O(log n + k) |

Where k is the number of results returned.

### Usage

B-trees are primarily used for:
- **Database systems**: Minimizing disk I/O
- **File systems**: Efficient file indexing
- **Large branching factor**: When m is large (e.g., m=256)

---

## Implementation Details

### AVL Tree Balancing

All statistic trees use AVL balancing to maintain O(log n) height.

#### Balance Factor

```cpp
template<typename T, typename V, typename S>
int balance(statistic_node<T, V, S>* tree) {
    return height(tree->left) - height(tree->right);
}
```

**Invariant**: Balance factor ∈ {-1, 0, 1} for all nodes.

#### Rotations

```cpp
// Right rotation
template<typename T, typename V, typename S>
statistic_node<T,V,S>* rebalance_right(statistic_node<T,V,S>* x) {
    auto y = x->left;
    auto B = y->right;

    y->right = x;
    x->left = B;

    // Update parent pointers and statistics
    x->update();
    y->update();

    return y;
}

// Left rotation (symmetric)
template<typename T, typename V, typename S>
statistic_node<T,V,S>* rebalance_left(statistic_node<T,V,S>* x);
```

#### Rebalancing After Modification

```cpp
template<typename T, typename V, typename S>
statistic_node<T,V,S>* rebalance(statistic_node<T,V,S>* x) {
    if(balance(x) < -1) {
        // Right-heavy
        if(balance(x->right) == 1)
            rebalance_right(x->right);  // Right-Left case
        x = rebalance_left(x);
    }
    else if(balance(x) > 1) {
        // Left-heavy
        if(balance(x->left) == -1)
            rebalance_left(x->left);  // Left-Right case
        x = rebalance_right(x);
    }

    x->update();  // Update statistics
    return rebalance(x->parent);  // Propagate upward
}
```

### Insertion Algorithm

```cpp
template<typename T, typename V, typename S>
statistic_node<T,V,S>* insert(
    statistic_node<T,V,S>* tree,
    const T& v,
    const V& data,
    bool or_assign = false
) {
    if(!tree) {
        return new statistic_node<T,V,S>(v, data);
    }

    auto p = lower_bound(tree, v);

    if(p == nullptr) {
        // Insert at rightmost position
        p = tree;
        while(p->right) p = p->right;
    }
    else if(or_assign && p->v == v) {
        // Update existing
        p->data = data;
        p->update();
        return rebalance(p);
    }
    else if(p->left) {
        // Find predecessor
        p = p->left;
        while(p->right) p = p->right;
    }

    // Create new node
    auto u = new statistic_node<T,V,S>(v, data, p);

    if(v <= p->v)
        p->left = u;
    else
        p->right = u;

    p->update();
    return rebalance(p);
}
```

**Time Complexity**: O(log n)

### Deletion Algorithm

```cpp
template<typename T, typename V, typename S>
std::pair<statistic_node<T,V,S>*, statistic_node<T,V,S>*>
extract(statistic_node<T,V,S>* tree, const T& v) {
    auto p = lower_bound(tree, v);

    if(!p->left) {
        // No left child: replace with right child
        auto w = p->parent ? p->parent : p->right;
        if(p->parent) {
            if(p->parent->left == p)
                p->parent->left = p->right;
            else
                p->parent->right = p->right;
        }
        if(p->right) p->right->parent = p->parent;

        p->right = nullptr;
        p->parent = nullptr;
        return {p, rebalance(w)};
    }
    else if(!p->left->right) {
        // Left child has no right: promote left child
        auto w = p->left;
        if(p->parent) {
            if(p->parent->left == p)
                p->parent->left = p->left;
            else
                p->parent->right = p->left;
        }
        if(p->right) p->right->parent = w;
        w->right = p->right;
        w->parent = p->parent;

        p->right = nullptr;
        p->left = nullptr;
        p->parent = nullptr;
        return {p, rebalance(w)};
    }
    else {
        // Replace with in-order predecessor
        auto u = p->left;
        while(u->right) u = u->right;

        auto s = u->parent;
        s->right = u->left;
        if(u->left) u->left->parent = s;

        std::swap(u->v, p->v);
        std::swap(u->data, p->data);

        u->left = nullptr;
        u->right = nullptr;
        u->parent = nullptr;
        return {u, rebalance(s)};
    }
}
```

**Time Complexity**: O(log n)

### Split and Merge

#### Split

```cpp
template<typename T, typename V, typename S>
std::pair<statistic_node<T,V,S>*, statistic_node<T,V,S>*>
split(statistic_node<T,V,S>* node, T threshold) {
    if(!node)
        return {nullptr, nullptr};

    if(node->v < threshold) {
        auto [L, R] = split(node->right, threshold);
        return {merge_with_root(node, node->left, L), R};
    }
    else {
        auto [L, R] = split(node->left, threshold);
        return {L, merge_with_root(node, R, node->right)};
    }
}
```

**Time Complexity**: O(log n)

#### Merge

```cpp
template<typename T, typename V, typename S>
statistic_node<T,V,S>* merge(
    statistic_node<T,V,S>* left,
    statistic_node<T,V,S>* right
) {
    if(!left) return right;

    // Extract largest from left
    statistic_node<T,V,S>* last = left;
    while(last->right) last = last->right;

    auto [root, L] = extract(left, last->v);
    return merge_with_root(root, L, right);
}
```

**Time Complexity**: O(log n)

---

## API Reference

### Common Operations (All Statistic Trees)

```cpp
namespace cp::data_structures::stats_trees {

// Insert key-value pair
template<typename T, typename V, typename S>
statistic_node<T,V,S>* insert(
    statistic_node<T,V,S>* tree,
    const T& v,
    const V& data,
    bool or_assign = false
);

// Insert or update
template<typename T, typename V, typename S>
statistic_node<T,V,S>* insert_or_assign(
    statistic_node<T,V,S>* tree,
    const T& v,
    const V& data
);

// Find node by key
template<typename T, typename V, typename S>
statistic_node<T,V,S>* find(
    statistic_node<T,V,S>* node,
    const T& v
);

// Get value by key (throws if not found)
template<typename T, typename V, typename S>
V& value_at(
    statistic_node<T,V,S>* node,
    const T& v
);

// Erase node by key
template<typename T, typename V, typename S>
statistic_node<T,V,S>* erase(
    statistic_node<T,V,S>* tree,
    const T& v
);

// Update value
template<typename T, typename V, typename S>
statistic_node<T,V,S>* update(
    statistic_node<T,V,S>* tree,
    const T& v,
    const V& data
);

// Destroy tree (free memory)
template<typename T, typename V, typename S>
void destroy(statistic_node<T,V,S>* node);

// Lower bound (>= v)
template<typename T, typename V, typename S>
statistic_node<T,V,S>* lower_bound(
    statistic_node<T,V,S>* tree,
    const T& v
);

// Upper bound (> v)
template<typename T, typename V, typename S>
statistic_node<T,V,S>* upper_bound(
    statistic_node<T,V,S>* tree,
    const T& v
);

// Reverse lower bound (<= v)
template<typename T, typename V, typename S>
statistic_node<T,V,S>* reverse_lower_bound(
    statistic_node<T,V,S>* tree,
    const T& v
);

// Reverse upper bound (< v)
template<typename T, typename V, typename S>
statistic_node<T,V,S>* reverse_upper_bound(
    statistic_node<T,V,S>* tree,
    const T& v
);

// Get smallest element
template<typename T, typename V, typename S>
statistic_node<T,V,S>* begin(
    statistic_node<T,V,S>* tree
);

// Get successor
template<typename T, typename V, typename S>
statistic_node<T,V,S>* next(
    statistic_node<T,V,S>* tree
);

// Get predecessor
template<typename T, typename V, typename S>
statistic_node<T,V,S>* prev(
    statistic_node<T,V,S>* tree
);

// Split tree at threshold
template<typename T, typename V, typename S>
std::pair<statistic_node<T,V,S>*, statistic_node<T,V,S>*>
split(
    statistic_node<T,V,S>* node,
    T threshold
);

// Merge two trees
template<typename T, typename V, typename S>
statistic_node<T,V,S>* merge(
    statistic_node<T,V,S>* left,
    statistic_node<T,V,S>* right
);

}
```

---

## Comparison with Other Structures

### vs. std::set / std::map

| Feature | std::set/map | Order Statistic Tree |
|---------|--------------|---------------------|
| Insert/Delete | O(log n) | O(log n) |
| Find | O(log n) | O(log n) |
| **Select kth** | O(n) | **O(log n)** |
| **Rank** | O(n) | **O(log n)** |
| **Range Sum** | O(n) | **O(log n)** |
| **Median** | O(n) | **O(log n)** |
| **Split/Merge** | O(n) | **O(log n)** |
| Memory | O(n) | O(n) |
| Implementation | Red-Black Tree | AVL Tree + Statistics |

**When to use Order Statistic Tree:**
- Need kth element queries
- Need rank queries
- Need range aggregates
- Need median efficiently
- Need split/merge operations

**When to use std::set/map:**
- Only need insert/delete/find
- Don't need order statistics
- Want standard library guarantees

### vs. Segment Tree

| Feature | Segment Tree | Statistic Tree |
|---------|--------------|----------------|
| Range Query | O(log n) | O(log n) |
| Point Update | O(log n) | O(log n) |
| **Insert/Delete** | Not supported | **O(log n)** |
| **Order Statistics** | Not supported | **O(log n)** |
| **Dynamic Size** | Fixed | **Dynamic** |
| Memory | O(n) | O(n) |

**When to use Statistic Tree:**
- Need dynamic insertions/deletions
- Need order statistics
- Tree size changes frequently

**When to use Segment Tree:**
- Fixed array size
- Only need range queries
- Don't need insertions/deletions

---

## Advanced Operations

### Custom Statistics

You can define custom statistics by implementing the statistic concept:

```cpp
// Custom: Track max value in subtree
struct max_value_stats {
    int size;
    int max_val;

    template<typename T>
    max_value_stats(const T& key, int value)
        : size(1), max_val(value) {}

    template<typename T>
    static void update(statistic_node<T, int, max_value_stats>* node) {
        node->statistic.size =
            (node->left ? node->left->statistic.size : 0) +
            1 +
            (node->right ? node->right->statistic.size : 0);

        node->statistic.max_val = node->data;
        if(node->left)
            node->statistic.max_val = std::max(
                node->statistic.max_val,
                node->left->statistic.max_val
            );
        if(node->right)
            node->statistic.max_val = std::max(
                node->statistic.max_val,
                node->right->statistic.max_val
            );
    }
};
```

### Persistent Statistic Trees

For persistent (immutable) versions:
- Use path copying on modifications
- Share unmodified subtrees between versions
- O(log n) space per update
- Access to all historical versions

### Range Update Operations

For lazy propagation:
- Add `lazy` field to statistics
- Propagate updates on-demand
- Supports range updates in O(log n)

---

## Best Practices

### 1. Choose the Right Statistic

```cpp
// Need order statistics only?
order_node<int>* tree;

// Need key sums?
key_sum_node_t<int, plus_t>* tree;

// Need value sums?
sum_node_t<int, double, plus_t>* tree;
```

### 2. Memory Management

```cpp
// Always destroy when done
destroy(tree);
tree = nullptr;

// Or use RAII wrapper
class TreeRAII {
    order_node<int>* tree;
public:
    TreeRAII() : tree(nullptr) {}
    ~TreeRAII() { destroy(tree); }
    // ...
};
```

### 3. Handle Empty Trees

```cpp
// Check before operations
if(tree) {
    auto median_node = median(tree);
}

// Or use size check
if(size(tree) > 0) {
    // Safe to operate
}
```

### 4. Use insert_or_assign for Updates

```cpp
// Update if exists, insert if not
tree = insert_or_assign(tree, key, new_value);

// Better than:
// if(find(tree, key))
//     tree = update(tree, key, new_value);
// else
//     tree = insert(tree, key, new_value);
```

---

## Performance Considerations

### Cache Efficiency

Statistic trees have:
- **Good locality**: Parent pointers enable upward traversal
- **Moderate memory**: O(n) nodes with small overhead
- **Pointer chasing**: May have cache misses on deep trees

### Comparison with Arrays

For static data:
- **Arrays**: O(1) access, O(n) insert/delete
- **Statistic trees**: O(log n) for all operations

Trade-off: Use arrays for static data, trees for dynamic.

### Optimization Tips

```cpp
// 1. Reserve memory (not applicable, uses new)

// 2. Batch operations
for(auto x : elements) {
    tree = insert(tree, x, {});
}
// Better than individual inserts with queries

// 3. Use appropriate value type
order_node<int, std::monostate>* tree;  // No value overhead

// 4. Profile before optimizing
// Measure actual performance needs
```

---

## Conclusion

CPLibrary's statistic trees provide powerful augmented BST functionality:

- **Order Statistic Trees**: Efficient kth element, rank, median queries
- **Key Sum Trees**: Range aggregates over keys
- **Value Sum Trees**: Range aggregates over values
- **B-Trees**: Multi-way trees for database systems

All with O(log n) operations and flexible generic design.

### Quick Selection Guide

**Use Order Statistic Tree when:**
- Need kth smallest/largest
- Need rank queries
- Need median
- Need dynamic insertions/deletions

**Use Key Sum Tree when:**
- Need range sums/products of keys
- Keys represent quantities to aggregate

**Use Value Sum Tree when:**
- Need range sums/products of values
- Values represent weights/counts

**Use B-Tree when:**
- Large branching factor needed
- Optimizing for disk I/O
- Database/file system applications

For more examples and integration, see:
- [EXAMPLES.md](../EXAMPLES.md) - Usage examples
- [API.md](../API.md) - Complete API reference
- [RANGE_QUERY_DATA_STRUCTURES.md](RANGE_QUERY_DATA_STRUCTURES.md) - Related structures

---

**Last Updated**: December 2025
**Documentation Author**: CPLibrary Documentation Team
**License**: See repository for license information
