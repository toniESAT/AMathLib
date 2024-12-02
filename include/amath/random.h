#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <math.h>
#include <random>
#include "utils.h"

namespace amath {

/**
 * @brief Uniform random integer generator
 * 
 * Generates random integers uniformly distributed in a specified range
 * using the Mersenne Twister engine (mt19937).
 */
struct RandomIntUniform {
   /**
    * @brief Constructs random generator with specified range
    * @param low Lower bound (inclusive)
    * @param high Upper bound (inclusive)
    */
   RandomIntUniform(int low, int high) : low(low), high(high) {
      gen = std::mt19937(rd());
      dist = std::uniform_int_distribution<int>(low, high);
   }

   /**
    * @brief Sets a specific seed for reproducibility
    * @param seed Seed value
    */
   void seed(unsigned int seed) { gen.seed(seed); };

   /**
    * @brief Changes the range of generated numbers
    * @param low New lower bound (inclusive)
    * @param high New upper bound (inclusive)
    */
   void setRange(int low, int high) { dist = std::uniform_int_distribution<int>(low, high); }

   /**
    * @brief Generates a random integer in the specified range
    * @return Random integer between low and high (inclusive)
    */
   int generate() { return dist(gen); }

   private:
   int low, high;                                   ///< Range bounds
   std::random_device rd;                          ///< Random device for seeding
   std::mt19937 gen;                              ///< Mersenne Twister generator
   std::uniform_int_distribution<int> dist;       ///< Integer distribution
};

/**
 * @brief Uniform random floating-point generator
 * 
 * Generates random floating-point numbers uniformly distributed in a specified range
 * using the Mersenne Twister engine (mt19937). Uses the scalar type defined in utils.h.
 */
struct RandomFloatUniform {
   /**
    * @brief Constructs random generator with specified range
    * @param low Lower bound (inclusive)
    * @param high Upper bound (inclusive)
    */
   RandomFloatUniform(scalar low, scalar high) : low(low), high(high) {
      gen = std::mt19937(rd());
      dist = std::uniform_real_distribution<scalar>(low, high);
   }

   /**
    * @brief Sets a specific seed for reproducibility
    * @param seed Seed value
    */
   void seed(unsigned int seed) { gen.seed(seed); };

   /**
    * @brief Changes the range of generated numbers
    * @param low New lower bound (inclusive)
    * @param high New upper bound (inclusive)
    */
   void setRange(scalar low, scalar high) {
      dist = std::uniform_real_distribution<scalar>(low, high);
   }

   /**
    * @brief Generates a random floating-point number in the specified range
    * @return Random number between low and high (inclusive)
    */
   scalar generate() { return dist(gen); }

   private:
   scalar low, high;                                    ///< Range bounds
   std::random_device rd;                              ///< Random device for seeding
   std::mt19937 gen;                                   ///< Mersenne Twister generator
   std::uniform_real_distribution<scalar> dist;        ///< Floating-point distribution
};

}  // namespace amath

#endif /* __RANDOM_H__ */