#ifndef __AMATH_BIT_H__
#define __AMATH_BIT_H__

#include <stdint.h>

namespace amath {

/**
 * @brief Checks if a number is a power of two
 *
 * Uses a bit manipulation trick to efficiently determine if a number is a power of two.
 * A power of two has exactly one bit set in its binary representation, and subtracting
 * one from it will set all lower bits to 1. The bitwise AND of these will be 0 only
 * for powers of two.
 *
 * Source: http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
 *
 * @param x Number to check
 * @return true if x is a power of two, false otherwise
 * @note Also returns false for x = 0
 */
bool IsPowerOfTwo(uint64_t x) { return (x != 0) && ((x & (x - 1)) == 0); }

}  // namespace amath

#endif /* __AMATH_BIT_H__ */