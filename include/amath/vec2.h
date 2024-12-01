#ifndef __VEC2_H__
#define __VEC2_H__

#include <math.h>
#include <stdio.h>

#include "common_defs.h"
#include "utils.h"

namespace amath {

struct Vec2 {
   scalar d[2];

   /*****************************
    ****  Vec2 constructors  ****
    *****************************/

   Vec2(scalar v0, scalar v1) : d{v0, v1} {};
   Vec2(scalar v) : d{v, v} {};
   Vec2() : Vec2(0, 0) {};

   static Vec2 nan() { return {nanf(""), nanf("")}; }

   /*****************************
    *****  Getters/Setters  *****
    *****************************/

   scalar x() const { return d[0]; }
   scalar y() const { return d[1]; }
   scalar &x() { return d[0]; }
   scalar &y() { return d[1]; }

   /*****************************
    *********  Methods  *********
    *****************************/

   scalar lengthSquared() const { return x() * x() + y() * y(); }
   scalar length() const { return sqrtf(lengthSquared()); }

   Vec2 normalized(scalar tolerance = AM_EPSILON) const {
      scalar l = length();
      if (!isAlmostZero(l, tolerance)) return {x() / l, y() / l};
      else return {0, 0};
   }

   bool isNormalized(const scalar tolerance = AM_EPSILON) const {
      return isAlmostZero(lengthSquared() - 1, tolerance);
   }

   Vec2 perpendicular() const { return {-y(), x()}; }

   inline scalar dot(const Vec2 &v) const { return x() * v.x() + y() * v.y(); }

   void print() const { printf("Vec2 [%.4f, %.4f]\n", x(), y()); }

   /*****************************
    *  Vec2 operator overloads  *
    *****************************/

   // Negate
   Vec2 operator-() const { return {-d[0], -d[1]}; }

   // Operations with scalars
   Vec2 operator+(const scalar k) const { return {x() + k, y() + k}; }
   Vec2 operator-(const scalar k) const { return {x() - k, y() - k}; }
   Vec2 operator*(const scalar k) const { return {x() * k, y() * k}; }
   Vec2 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      return *this;
   }
   Vec2 operator/(const scalar k) const { return {x() / k, y() / k}; }
   Vec2 &operator/=(const scalar k) {
      d[0] /= k;
      d[1] /= k;
      return *this;
   }

   // Operations with vectors
   Vec2 operator+(const Vec2 &v) const { return {x() + v.x(), y() + v.y()}; }
   Vec2 operator-(const Vec2 &v) const { return {x() - v.x(), y() - v.y()}; }
   Vec2 &operator+=(Vec2 &v) {
      x() += v.x();
      y() += v.y();
      return *this;
   }
   Vec2 &operator-=(Vec2 &v) {
      x() -= v.x();
      y() -= v.y();
      return *this;
   }

   // Comparison
   bool operator==(const Vec2 &v) const {
      return fabs(x() - v.x()) < AM_EPSILON && fabs(y() - v.y()) < AM_EPSILON;
   }
   bool operator!=(const Vec2 &v) const {
      return fabs(x() - v.x()) > AM_EPSILON || fabs(y() - v.y()) > AM_EPSILON;
   }

   // Access
   scalar operator[](const size_t i) const { return d[i]; }
   scalar &operator[](const size_t i) { return d[i]; }
};

}  // namespace amath

#endif /* __VEC2_H__ */
