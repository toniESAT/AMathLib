#ifndef __COMMON_DEFS_H__
#define __COMMON_DEFS_H__

#include <float.h>
#include <stdint.h>

namespace amath {

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t u8;
typedef int64_t s64;
typedef int32_t s32;
typedef int8_t s8;

#ifdef AM_DOUBLE_PRECISON
typedef double scalar;
#define AM_EPSILON DBL_EPSILON
#else
typedef float scalar;
#define AM_EPSILON FLT_EPSILON
#endif

const scalar kSmallNumber = (scalar)1.e-8;
const scalar PI = (scalar)3.141592653589793;

}  // namespace amath

#endif /* __COMMON_DEFS_H__ */
