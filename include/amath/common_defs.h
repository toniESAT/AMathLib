#ifndef __COMMON_DEFS_H__
#define __COMMON_DEFS_H__

#include <float.h>
#include <stdint.h>

#define u64 uint64_t
#define u32 uint32_t
#define u8 uint8_t
#define s64 int64_t
#define s32 int32_t
#define s8 int8_t

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

namespace amath {

#ifdef AM_DOUBLE_PRECISON
typedef double scalar;
#define AM_EPSILON DBL_EPSILON
#else
typedef float scalar;
#define AM_EPSILON FLT_EPSILON
#endif

const scalar kSmallNumber = 1.e-8f;
const scalar PI = 3.141592653589793;

}  // namespace amath

#endif /* __COMMON_DEFS_H__ */
