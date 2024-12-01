#ifndef __BIT_H__
#define __BIT_H__

#include <stdint.h>

namespace amath {
// http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
bool IsPowerOfTwo(uint64_t x) { return (x != 0) && ((x & (x - 1)) == 0); }
}  // namespace amath

#endif /* __BIT_H__ */
