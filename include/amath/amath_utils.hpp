#pragma once

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
constexpr float kSmallNumber = 1.e-8f;

constexpr float PI = 3.141592653589793;
constexpr float Rad2Deg(float a) { return a * 180.f / PI; }
constexpr float Deg2Rad(float a) { return a * PI / 180.f; }

bool is_almost_zero(float x, float tolerance = FLT_EPSILON) { return fabsf(x) < tolerance; }
bool is_almost_zero(double x, double tolerance = DBL_EPSILON) { return fabs(x) < tolerance; }

template <typename T>
std::vector<size_t> argsort(const std::vector<T> &array, const bool descending = false) {
   std::vector<size_t> indices(array.size());
   for (size_t i = 0; i < array.size(); i++) indices[i] = i;
   std::sort(indices.begin(), indices.end(), [&array, descending](size_t left, size_t right) -> bool {
      if (!descending) return array[left] < array[right];
      else return array[left] > array[right];
   });

   return indices;
}

} // namespace amath