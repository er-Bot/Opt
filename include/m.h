#pragma once

#include <cmath>
#include <cstdlib>

/***             Math Constants               ***/
#define GOLDEN_RATIO 1.6180339887499
#define EPSILON 1e-3

#define INFTY INT64_MAX
#define FROBINUS INT64_MAX-1

/***             Basic Functions               ***/
#define _min(__a, __b)      (__a > __b) ? __b : __a
#define _max(__a, __b)      (__a < __b) ? __b : __a
#define _abs(__a)           (__a < 0) ? -__a : __a
#define _sgn(__a)           (__a < 0) ? -1 : 1
#define _swap(__a, __b)     {auto __t = __b; __b = __a; __a = __t;} 
#define _rnd()             ((double) rand() / (RAND_MAX - 1))
#define _rand(__a, __b)     _rnd() * (__b - __a) + __a

size_t fib(size_t n);