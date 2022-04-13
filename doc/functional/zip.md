# zip

```c++
#include "function/zip.h"
```



## 1. Introduction

`zip` is a well known and powerful function in Python that from an iterators $I_1,\dots,I_n$ over iterables $A_1,\dots,A_n$ respectively (that have the same size), it can build an iterator $\bold{I}$ over $\bold{A}=(A_1,\dots,A_n)$.

## 2. Function

`zip`'s iterator can function both by (possibly const r/l value) reference, or by copy. This can be controlled by the functions:

- `zip` for default r-value reference
- `zip_copy` for copying values
- `zip_const` for const l-value reference