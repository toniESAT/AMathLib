#ifndef __UTILS_H__
#define __UTILS_H__

#include <float.h>
#include <vector>
#include <algorithm>
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

inline scalar Rad2Deg(scalar a) { return a * 180.f / PI; }
inline scalar Deg2Rad(scalar a) { return a * PI / 180.f; }

inline bool isAlmostZero(scalar x, scalar tolerance = FLT_EPSILON) { return fabsf(x) < tolerance; }

template <typename T>
inline std::vector<size_t> argsort(const std::vector<T> &array, const bool descending = false) {
   std::vector<size_t> indices(array.size());
   for (size_t i = 0; i < array.size(); i++) indices[i] = i;
   std::sort(
       indices.begin(), indices.end(), [&array, descending](size_t left, size_t right) -> bool {
          if (!descending) return array[left] < array[right];
          else return array[left] > array[right];
       });

   return indices;
}

}  // namespace amath

#endif /* __UTILS_H__ */
