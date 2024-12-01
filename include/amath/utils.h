#ifndef __UTILS_H__
#define __UTILS_H__

#include <float.h>
#include <stdint.h>

namespace amath {

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t u8;
typedef int64_t s64;
typedef int32_t s32;
typedef int8_t s8;

#ifdef AM_DOUBLE_PRECISONs
typedef double scalar;
#define AM_EPSILON DBL_EPSILON
#else
typedef float scalar;
#define AM_EPSILON FLT_EPSILON
#endif

const scalar kSmallNumber = (scalar)1.e-8;
const scalar PI = (scalar)3.141592653589793;

inline scalar Rad2Deg(const scalar a) { return a * (scalar)180.0 / PI; }
inline scalar Deg2Rad(const scalar a) { return a * PI / (scalar)180.0; }

inline bool almostZero(const scalar x, const scalar tolerance = AM_EPSILON) {
   return fabs(x) < AM_EPSILON;
}
inline bool almostEqual(const scalar x, const scalar y, const scalar tolerance = AM_EPSILON) {
   return fabs(x - y)  < AM_EPSILON;
}

}  // namespace amath

#endif /* __UTILS_H__ */
