#pragma once

#include <cfloat>

namespace amath {
constexpr float kSmallNumber = 1.e-8f;

constexpr float PI = 3.14159265359;

bool is_almost_zero(float x, float tolerance = FLT_EPSILON) { return fabsf(x) < tolerance; }
bool is_almost_zero(double x, double tolerance = DBL_EPSILON) { return fabs(x) < tolerance; }

} // namespace amath