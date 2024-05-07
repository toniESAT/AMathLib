#pragma once

#include <cfloat>
#include <vector>
#include <algorithm>

namespace amath {
constexpr float kSmallNumber = 1.e-8f;

constexpr float PI = 3.14159265359;

bool is_almost_zero(float x, float tolerance = FLT_EPSILON) { return fabsf(x) < tolerance; }
bool is_almost_zero(double x, double tolerance = DBL_EPSILON) { return fabs(x) < tolerance; }

template <typename T>
std::vector<size_t> argsort(const std::vector<T> &array, const bool descending = false) {
   std::vector<size_t> indices(array.size());
   for (size_t i = 0; i < array.size(); i++) indices[i] = i;
   std::sort(indices.begin(), indices.end(), [&array, descending](int left, int right) -> bool {
      if (!descending) return array[left] < array[right];
      else return array[left] > array[right];
   });

   return indices;
}

} // namespace amath