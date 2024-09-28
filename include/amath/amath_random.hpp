#pragma once

#include <math.h>
#include <random>
#include <immintrin.h>
#define FMT_HEADER_ONLY
#include "fmt/format.h"

namespace amath {

struct RandomIntUniform {
   RandomIntUniform(int low, int high) : low(low), high(high) {
      gen = std::mt19937(rd());
      dist = std::uniform_int_distribution<int>(low, high);
   }

   void seed(unsigned int seed) { gen.seed(seed); };

   void set_range(int low, int high) { dist = std::uniform_int_distribution<int>(low, high); }

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

   void set_range(float low, float high) {
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