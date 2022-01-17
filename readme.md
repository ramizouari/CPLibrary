# Competitive Programming Library
## 1. Rationale
During my Competitive Programming's journey, I spent a good time exploring novel approaches on many interesting problems.
This lead me to learn many new algorithms and data structures, and even incorporate my knowledge from different fields (Linear Algebra, Abstract Algebra...) to 
solve many hard problems. So I decided to write them in this repository.

I chose to write that acquired knowledge on this repository so that I can share my wonderful journey with any interested person.



## 2. Content
This repository contains a list of non-trivial algorithms & data structures that I have used at least once in competitive programming settings.
It is written with emphasis on readability, zero-overhead, and also correctness and mathematical rigour.

The following list of topics are present on this library:
1. Data Structures: 
	- Segment Tree over an associative binary operation
	- Fenwick Tree over a group
	- Sparse Table over an assoviative idempotent binary operation
	- Order Statistic Tree
	- Statistic Trees with Statistic acting on keys and/or values (over an associative binary operation)
2. String (To Be Completed):
	- Rabin-Karp probabilistic string matching algorithm.
	- ...
3. Graph (To Be Done)
4. Abstract Algebra:
	- Generic definition of usual binary operations
	- Fast Exponentiation over a monoid <img src="https://render.githubusercontent.com/render/math?math=\mathcal{M}">
	- Extended Euclidean Algorithm, and Bezout Coefficient over an integral domain <img src="https://render.githubusercontent.com/render/math?math=\mathcal{I}">
	- Rational Extension: field of rationals over an integral domain <img src="https://render.githubusercontent.com/render/math?math=\mathcal{I}">
	- dynamic/static ring extension of a commutative ring
	- Quadratic extension of a commutative ring
	- Conjugate element of a quadratic extension
	- Consistent division between two elements of a ring extension if the base ring is an integral domain
	- Inverse of an element of a ring extension if the base ring is a field
5. Modular Arithmetic
	- Support for static and dynamic cyclic elements
	- Modular Arithmetic
	- Primitive Roots of unity
	- Discrete logarithm
	- Modular square root
	- Legendre Symbol
	- Linear time inverse table: <img src="https://render.githubusercontent.com/render/math?math=1^{-1},\dots,n^{-1}"> in 
	<img src="https://render.githubusercontent.com/render/math?math=\mathcal{O}(n)">.
6. Linear Algebra:
	- dynamic and static vectors over commutative rings
	- dynamic and static matrices over commutative rings
	- Matrix arithmetic
	- Solving of a linear system over a field
	- Calculating matrix determinant
	- Calculating the characteristic polynomial of a matrix over a field or an integral domain.
7. Polynomials(To Be Complete):
	- Polynomials and sparse polynomials over a commutative ring
	- Polynomials arithmetic
	- subquadratic polynomial multiplication for commutative rings
	- Newton Interpolation
8. FFT(To Be Complete):
	- FFT over <img src="https://render.githubusercontent.com/render/math?math=\mathbb{C}">
	- FFT over <img src="https://render.githubusercontent.com/render/math?math=\mathbb{F}_p"> 
	for vectors with size <img src="https://render.githubusercontent.com/render/math?math=n | p">
	- FFT2 over <img src="https://render.githubusercontent.com/render/math?math=\mathbb{C}">
	- FFT2 over <img src="https://render.githubusercontent.com/render/math?math=\mathbb{F}_p"> for vectors with size 
	<img src="https://render.githubusercontent.com/render/math?math=2^k | p">
	- Fast Hadamard Transform over a ring <img src="https://render.githubusercontent.com/render/math?math=\mathcal{R}">
	- Multidimensional FFT over <img src="https://render.githubusercontent.com/render/math?math=\mathbb{C}">
	- Multidimensional FFT over <img src="https://render.githubusercontent.com/render/math?math=\mathbb{F}_p"> 
	if the tensor is of shape <img src="https://render.githubusercontent.com/render/math?math=(n_1,\dots,n_k)"> and 
	<img src="https://render.githubusercontent.com/render/math?math=n_i | p \quad \forall i\in\{1,\dots,k\}"> 
	- Fast polynomial multiplication over <img src="https://render.githubusercontent.com/render/math?math=\mathbb{C},\mathbb{R},\mathbb{Z},\mathbb{F}_p">
9. Fast Polynomials(To Be Done)
10. Number Theory:
	- Prime Numbers sieve in <img src="https://render.githubusercontent.com/render/math?math=\mathcal{O}(n\ln n)">
	- Implementation of the most used multiplicative functions
	- Chinese Remainder Theorem
	- Rabin-Miller primality test
	- Rho algorithm for factorisation
11.	Optimization:
	- Simplex Method

## 3. Additional Content
The following list of topics was not intended for CP, but I chose to include it as a demonstration to the power of this library:
1.	Analysis:
	- First <img src="https://render.githubusercontent.com/render/math?math=n"> 
	terms of <img src="https://render.githubusercontent.com/render/math?math=\exp,\log"> over a suitable commutative ring
	- <img src="https://render.githubusercontent.com/render/math?math=\exp,\log"> over a Banach algebra $\mathcal{A}$
	- Newton-Raphson Method over a <img src="https://render.githubusercontent.com/render/math?math=\mathbb{K}">-vector space, 
	with <img src="https://render.githubusercontent.com/render/math?math=\mathbb{K}\in\{\mathbb{R},\mathbb{C}\}">.
2.	Topology:
	- Abstract metrics, norms and inner products
	- Generic definition for the usual metrics, norms and inner products.
	- Derivation of any multivariable function between two <img src="https://render.githubusercontent.com/render/math?math=\mathbb{K}">-vector 
	spaces, with <img src="https://render.githubusercontent.com/render/math?math=\mathbb{K}\in\{\mathbb{R},\mathbb{C}\}">.
3.	Functional:
	- pointwise unary operators: calculate a scalar function pointwise on all elements of a matrix, vector or even an iterable
	- pointwise binary operators
	- foreach over iterables
	- pointwise reduce: reduce all elements under an associative binary operator
	- pointwise aggregate: transform each element of a matrix/vector to a given element, and reduce the resulting matrix/vector under
		an associative operator
4.	Generic:
	- zip function: iterate simultaneously over iterables having the same size.
5.	Basic Machine Learning Support:
	- Linear Regression
	- Logistic Regression
	- Multinomial Logistic Regression
	- K-Neighbours Regression
	- K-Neighbours Classification
6.	Order:
	- Order closure of a totally ordered set <img src="https://render.githubusercontent.com/render/math?math=S">.
	- Algebraic operations on an order closure of a group/ring having a total order.
7. B-Tree:
	- B-Trees over an underlying Order Statistic Tree.
	- Guaranteed to have <img src="https://render.githubusercontent.com/render/math?math=\mathcal{O}(\ln n)"> performance, 
	even for large values of <img src="https://render.githubusercontent.com/render/math?math=m">.

## 4. Remarks
This is only a tiny subset of the computer science & mathematical litterature, so I recommend you to frequently read and discover new topics.