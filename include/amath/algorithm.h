#ifndef __ALGORITHM_H__
#define __ALGORITHM_H__

#include <vector>
#include <algorithm>

namespace amath {

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

#endif /* __ALGORITHM_H__ */