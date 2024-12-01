#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <math.h>
#include <random>

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
   RandomFloatUniform(float low, float high) : low(low), high(high) {
      gen = std::mt19937(rd());
      dist = std::uniform_real_distribution<float>(low, high);
   }

   void seed(unsigned int seed) { gen.seed(seed); };

   void setRange(float low, float high) {
      dist = std::uniform_real_distribution<float>(low, high);
   }

   float generate() { return dist(gen); }

   private:
   float low, high;
   std::random_device rd;
   std::mt19937 gen;
   std::uniform_real_distribution<float> dist;
};

}  // namespace amath

#endif /* __RANDOM_H__ */
