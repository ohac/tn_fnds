#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <xmmintrin.h>
#include <pmmintrin.h>
using std::min;
using std::max;
typedef uint32_t DWORD;
#define _aligned_malloc(a,b) aligned_alloc(b,a)
#define _aligned_free(a)
