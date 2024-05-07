#pragma once

#include <cmath>
#include <random>
#include <immintrin.h>
#define FMT_HEADER_ONLY
#include <fmt/format.h>

namespace amath {

struct RandomIntUniform {
   int low, high;
   std::random_device rd;
   std::mt19937 gen;
   std::uniform_int_distribution<int> dist;

   RandomIntUniform(int low, int high) : low(low), high(high) {
      gen = std::mt19937(rd());
      dist = std::uniform_int_distribution<int>(low, high);
   }

   void seed(unsigned int seed) { gen.seed(seed); };

   void set_range(int low, int high) { dist = std::uniform_int_distribution<int>(low, high); }

   int generate() { return dist(gen); }
};

struct RandomFloatUniform {
   float low, high;
   std::random_device rd;
   std::mt19937 gen;
   std::uniform_real_distribution<float> dist;

   RandomFloatUniform(float low, float high) : low(low), high(high) {
      gen = std::mt19937(rd());
      dist = std::uniform_real_distribution<float>(low, high);
   }

   void seed(unsigned int seed) { gen.seed(seed); };

   void set_range(float low, float high) { dist = std::uniform_real_distribution<float>(low, high); }

   float generate() { return dist(gen); }
};

} // namespace amath