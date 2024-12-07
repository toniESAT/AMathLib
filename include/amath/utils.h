#ifndef __AMATH_UTILS_H__
#define __AMATH_UTILS_H__

#include <float.h>
#include <math.h>
#include <stdint.h>

namespace amath {

/**
 * @name Integer Type Aliases
 * Fixed-width integer type definitions for clarity and portability
 */
///@{
typedef uint64_t u64;  ///< 64-bit unsigned integer
typedef uint32_t u32;  ///< 32-bit unsigned integer
typedef uint8_t u8;    ///< 8-bit unsigned integer
typedef int64_t s64;   ///< 64-bit signed integer
typedef int32_t s32;   ///< 32-bit signed integer
typedef int8_t s8;     ///< 8-bit signed integer
///@}

/**
 * @brief Configurable floating-point precision type
 *
 * When AM_DOUBLE_PRECISION is defined, uses double precision;
 * otherwise uses single precision (float).
 */
#ifdef AM_DOUBLE_PRECISON
typedef double scalar;
#define AM_EPSILON DBL_EPSILON  ///< Machine epsilon for double precision
#else
typedef float scalar;
#define AM_EPSILON FLT_EPSILON  ///< Machine epsilon for single precision
#endif

/** @brief Small number constant for numerical comparisons */
const scalar kSmallNumber = (scalar)1.e-8;

/** @brief Mathematical constant π */
const scalar PI = (scalar)3.141592653589793;

/**
 * @brief Converts angle from radians to degrees
 * @param a Angle in radians
 * @return Angle in degrees
 */
inline scalar Rad2Deg(const scalar a) { return a * (scalar)180.0 / PI; }

/**
 * @brief Converts angle from degrees to radians
 * @param a Angle in degrees
 * @return Angle in radians
 */
inline scalar Deg2Rad(const scalar a) { return a * PI / (scalar)180.0; }

/**
 * @brief Checks if a value is approximately zero
 * @param x Value to check
 * @param tolerance Maximum difference from zero (defaults to machine epsilon)
 * @return true if |x| < tolerance
 */
inline bool almostZero(const scalar x, const scalar tolerance = AM_EPSILON) {
   return fabs(x) < AM_EPSILON;
}

/**
 * @brief Checks if two values are approximately equal
 * @param x First value
 * @param y Second value
 * @param tolerance Maximum allowed difference (defaults to machine epsilon)
 * @return true if |x - y| < tolerance
 */
inline bool almostEqual(const scalar x, const scalar y, const scalar tolerance = AM_EPSILON) {
   return fabs(x - y) < AM_EPSILON;
}

}  // namespace amath

#endif /* __AMATH_UTILS_H__ */