#ifndef __VEC3_H__
#define __VEC3_H__

#include <math.h>
#include <stdio.h>

#include "common_defs.h"
#include "utils.h"

namespace amath {

struct Vec3 {
   scalar d[3];

   /*****************************
    ****  Vec3 constructors  ****
    *****************************/

   Vec3(scalar v0, scalar v1, scalar v2) : d{v0, v1, v2} {};
   Vec3(scalar v) : d{v, v, v} {};
   Vec3() : Vec3(0, 0, 0) {};

   static Vec3 up() { return {0, 1, 0}; }
   static Vec3 nan() { return {nanf(""), nanf(""), nanf("")}; }

   /*****************************
    *****  Getters/Setters  *****
    *****************************/

   scalar x() const { return d[0]; }
   scalar y() const { return d[1]; }
   scalar z() const { return d[2]; }
   scalar &x() { return d[0]; }
   scalar &y() { return d[1]; }
   scalar &z() { return d[2]; }

   /*****************************
    *********  Methods  *********
    *****************************/

   scalar lengthSquared() const { return x() * x() + y() * y() + z() * z(); }
   scalar length() const { return sqrtf(lengthSquared()); }

   Vec3 normalized(scalar tolerance = AM_EPSILON) const {
      scalar l = length();
      if (!isAlmostZero(l, tolerance)) return {x() / l, y() / l, z() / l};
      else return {0, 0, 0};
   }

   bool isNormalized(const scalar tolerance = AM_EPSILON) const {
      return isAlmostZero(lengthSquared() - 1, tolerance);
   }

   scalar dot(const Vec3 &v) const { return x() * v.x() + y() * v.y() + z() * v.z(); }

   Vec3 cross(const Vec3 &v) const {
      return {y() * v.z() - z() * v.y(), z() * v.x() - x() * v.z(), x() * v.y() - y() * v.x()};
   }

   void print() const { printf("Vec3 [%.4f,%.4f, %.4f]\n", x(), y(), z()); }

   /*****************************
    *  Vec3 operator overloads  *
    *****************************/
   // Negate
   Vec3 operator-() const { return {-d[0], -d[1], -d[2]}; }

   // Operations with scalars
   Vec3 operator+(const scalar k) const { return {x() + k, y() + k, z() + k}; }
   Vec3 operator-(const scalar k) const { return {x() - k, y() - k, z() - k}; }
   Vec3 operator*(const scalar k) const { return {x() * k, y() * k, z() * k}; }
   Vec3 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      return *this;
   }
   Vec3 operator/(const scalar k) const { return {x() / k, y() / k, z() / k}; }
   Vec3 &operator/=(const scalar k) {
      d[0] /= k;
      d[1] /= k;
      d[2] /= k;
      return *this;
   }

   // Operations with vectors
   Vec3 operator+(const Vec3 &v) const { return {x() + v.x(), y() + v.y(), z() + v.z()}; }
   Vec3 operator-(const Vec3 &v) const { return {x() - v.x(), y() - v.y(), z() - v.z()}; }
   Vec3 &operator+=(Vec3 &v) {
      x() += v.x();
      y() += v.y();
      z() += v.y();
      return *this;
   }
   Vec3 &operator-=(Vec3 &v) {
      x() -= v.x();
      y() -= v.y();
      z() -= v.y();
      return *this;
   }

   // Comparison
   bool operator==(const Vec3 &v) const {
      return fabs(x() - v.x()) < AM_EPSILON && fabs(y() - v.y()) < AM_EPSILON &&
             fabs(z() - v.z()) < AM_EPSILON;
   }
   bool operator!=(const Vec3 &v) const {
      return fabs(x() - v.x()) > AM_EPSILON || fabs(y() - v.y()) > AM_EPSILON ||
             fabs(z() - v.z()) > AM_EPSILON;
   }

   // Access
   scalar operator[](const size_t i) const { return d[i]; }
   scalar &operator[](const size_t i) { return d[i]; }
};

}  // namespace amath

#endif /* __VEC3_H__ */
