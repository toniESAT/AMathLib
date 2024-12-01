#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <math.h>
#include <random>

#include "utils.h"

namespace amath {

struct RandomIntUniform {
   RandomIntUniform(int low, int high) : low(low), high(high) {
      gen = std::mt19937(rd());
      dist = std::uniform_int_distribution<int>(low, high);
   }

   void seed(unsigned int seed) { gen.seed(seed); };

   void setRange(int low, int high) { dist = std::uniform_int_distribution<int>(low, high); }

   int generate() { return dist(gen); }

   private:
   int low, high;
   std::random_device rd;
   std::mt19937 gen;
   std::uniform_int_distribution<int> dist;
};

struct RandomFloatUniform {
   RandomFloatUniform(scalar low, scalar high) : low(low), high(high) {
      gen = std::mt19937(rd());
      dist = std::uniform_real_distribution<scalar>(low, high);
   }

   void seed(unsigned int seed) { gen.seed(seed); };

   void setRange(scalar low, scalar high) {
      dist = std::uniform_real_distribution<scalar>(low, high);
   }

   scalar generate() { return dist(gen); }

   private:
   scalar low, high;
   std::random_device rd;
   std::mt19937 gen;
   std::uniform_real_distribution<scalar> dist;
};

}  // namespace amath

#endif /* __RANDOM_H__ */
