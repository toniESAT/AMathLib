#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <algorithm>

#include "common_defs.h"

namespace amath {

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
