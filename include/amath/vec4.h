#ifndef __VEC4_H__
#define __VEC4_H__

#include <math.h>
#include <stdio.h>

#include "utils.h"
#include "utils.h"

namespace amath {
struct Vec4 {
   scalar d[4];

   /*****************************
    ****  Vec4 constructors  ****
    *****************************/

   Vec4(scalar v0, scalar v1, scalar v2, scalar v3) : d{v0, v1, v2, v3} {};
   Vec4(scalar v) : d{v, v, v, v} {};
   Vec4() : Vec4(0, 0, 0, 0) {};

   static Vec4 up() { return {0, 1, 0, 0}; }
   static Vec4 nan() { return {nanf(""), nanf(""), nanf(""), nanf("")}; }

   /*****************************
    *****  Getters/Setters  *****
    *****************************/

   scalar x() const { return d[0]; }
   scalar y() const { return d[1]; }
   scalar z() const { return d[2]; }
   scalar w() const { return d[3]; }
   scalar &x() { return d[0]; }
   scalar &y() { return d[1]; }
   scalar &z() { return d[2]; }
   scalar &w() { return d[3]; }

   /*****************************
    *********  Methods  *********
    *****************************/

   scalar lengthSquared() const { return x() * x() + y() * y() + z() * z() + w() * w(); }
   scalar length() const { return sqrtf(lengthSquared()); }

   Vec4 normalized(scalar tolerance = AM_EPSILON) const {
      scalar k = 1 / length();
      // If homogeneous vector (w==0), normalize to length = 1
      if (almostZero(d[3]))
         //  if (almostZero(l, tolerance)) return {0, 0, 0, 0};  // TODO: Prevent div by 0 ???
         return {x() * k, y() * k, z() * k, 0};
      // If homogeneous point (w!=1), normalize to w = 1
      return *this * (1.f / d[3]);
   }

   bool isNormalized(const scalar tolerance = AM_EPSILON) {
      return almostZero(lengthSquared() - 1, tolerance);
   }

   scalar dot(const Vec4 &v) const { return x() * v.x() + y() * v.y() + z() * v.z() + w() * v.w(); }

   Vec4 cross(const Vec4 &v) const {
      return {y() * v.z() - z() * v.y(), z() * v.x() - x() * v.z(), x() * v.y() - y() * v.x(), 0};
   }

   void print() const { printf("Vec4 [%.4f, %.4f, %.4f, %.4f]\n", x(), y(), z(), w()); }

   /*****************************
    *  Vec4 operator overloads  *
    *****************************/

   // Negate
   Vec4 operator-() const { return {-d[0], -d[1], -d[2], -d[3]}; }

   // Operations with scalars
   Vec4 operator+(const scalar k) const { return {x() + k, y() + k, z() + k, w() + k}; }
   Vec4 operator-(const scalar k) const { return {x() - k, y() - k, z() - k, w() - k}; }
   Vec4 operator*(const scalar k) const { return {x() * k, y() * k, z() * k, w() * k}; }
   Vec4 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      return *this;
   }
   Vec4 operator/(const scalar k) const { return {x() / k, y() / k, z() / k, w() / k}; }
   Vec4 &operator/=(const scalar k) {
      d[0] /= k;
      d[1] /= k;
      d[2] /= k;
      d[3] /= k;
      return *this;
   }

   // Operations with vectors
   Vec4 operator-(const Vec4 &v) const {
      return {x() - v.x(), y() - v.y(), z() - v.z(), w() - v.w()};
   }
   Vec4 operator+(const Vec4 &v) const {
      return {x() + v.x(), y() + v.y(), z() + v.z(), w() + v.w()};
   }
   Vec4 &operator+=(Vec4 &v) {
      x() += v.x();
      y() += v.y();
      z() += v.z();
      w() += v.w();
      return *this;
   }
   Vec4 &operator-=(Vec4 &v) {
      x() -= v.x();
      y() -= v.y();
      z() -= v.z();
      w() -= v.w();
      return *this;
   }

   // Comparison
   bool operator==(const Vec4 &v) const {
      return fabs(x() - v.x()) < AM_EPSILON && fabs(y() - v.y()) < AM_EPSILON &&
             fabs(z() - v.z()) < AM_EPSILON && fabs(w() - v.w()) < AM_EPSILON;
   }
   bool operator!=(const Vec4 &v) const {
      return fabs(x() - v.x()) > AM_EPSILON || fabs(y() - v.y()) > AM_EPSILON ||
             fabs(z() - v.z()) > AM_EPSILON || fabs(w() - v.w()) > AM_EPSILON;
   }

   // Access
   scalar operator[](const size_t i) const { return d[i]; }
   scalar &operator[](const size_t i) { return d[i]; }
};

}  // namespace amath

#endif /* __VEC4_H__ */
