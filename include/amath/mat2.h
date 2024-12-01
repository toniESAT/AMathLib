#ifndef __MAT2_H__
#define __MAT2_H__

#include <math.h>
#include <stdio.h>
#include <vector>

#include "utils.h"

#include "vec2.h"

namespace amath {

struct Mat2 {
   scalar d[4];

   /*****************************
    ****  Mat2 constructors  ****
    *****************************/

   Mat2(scalar m0, scalar m1, scalar m2, scalar m3) : d{m0, m1, m2, m3} {};
   Mat2(scalar k) : Mat2(k, k, k, k) {}
   Mat2() : d{0} {};
   Mat2(const Vec2 &v0, const Vec2 &v1) {
      setCol(0, v0);
      setCol(1, v1);
   }

   static Mat2 identity() { return {1, 0, 0, 1}; }
   static Mat2 nan() { return Mat2(nanf("")); }

   /*****************************
    *****  Getters/Setters  *****
    *****************************/

   Vec2 getCol(const size_t j) const { return {d[2 * j], d[2 * j + 1]}; }
   Vec2 getRow(const size_t i) const { return {d[i], d[i + 2]}; }
   void setCol(const size_t j, const Vec2 &v) {
      d[2 * j] = v[0];
      d[2 * j + 1] = v[1];
   }
   void setRow(const size_t i, const Vec2 &v) {
      d[i] = v[0];
      d[i + 2] = v[1];
   }

   /*****************************
    *********  Methods  *********
    *****************************/

   static size_t size() { return 4; }

   scalar det() const { return d[0] * d[3] - d[1] * d[2]; }

   Mat2 transposed() const { return {d[0], d[2], d[1], d[3]}; };

   std::vector<Vec2> transformPoints(const std::vector<Vec2> &points) const {
      std::vector<Vec2> transformed_points;
      transformed_points.reserve(points.size());
      for (auto p : points) transformed_points.push_back(*this * p);
      return transformed_points;
   }

   void print() const { printf("Mat2:\n[%.4f][%.4f]\n[%.4f][%.4f]", d[0], d[2], d[1], d[3]); }

   /*****************************
    *  Mat2 operator overloads  *
    *****************************/

   // Access
   scalar operator[](size_t i) const { return d[i]; }
   scalar &operator[](size_t i) { return d[i]; }

   scalar operator()(size_t i, size_t j) const {
      if (i >= 0 && i < 2 && j >= 0 && j < 2) return d[i * 2 + j];
   }
   scalar &operator()(size_t i, size_t j) {
      if (i >= 0 && i < 2 && j >= 0 && j < 2) return d[i * 2 + j];
   }

   // Negate
   Mat2 operator-() const { return {-d[0], -d[1], -d[2], -d[3]}; }

   // Operations with scalars
   Mat2 operator+(const scalar k) const { return {d[0] + k, d[1] + k, d[2] + k, d[3] + k}; }
   Mat2 operator-(const scalar k) const { return *this + (-k); }

   Mat2 &operator+=(const scalar k) {
      d[0] += k;
      d[1] += k;
      d[2] += k;
      d[3] += k;

      return *this;
   }
   Mat2 &operator-=(const scalar k) {
      *this += (-k);
      return *this;
   }

   Mat2 operator*(const scalar k) const { return {d[0] * k, d[1] * k, d[2] * k, d[3] * k}; }
   Mat2 operator/(const scalar k) const { return *this * (1 / k); }

   Mat2 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      return *this;
   }
   Mat2 &operator/=(const scalar k) {
      *this *= (1 / k);
      return *this;
   }

   // Operations with vectors
   Vec2 operator*(const Vec2 &v) const {
      return {d[0] * v.x() + d[2] * v.y(), d[1] * v.x() + d[3] * v.y()};
   }

   // Operations with other matrices
   Mat2 operator+(const Mat2 &m) const {
      return {d[0] + m.d[0], d[1] + m.d[1], d[2] + m.d[2], d[3] + m.d[3]};
   }
   Mat2 operator-(const Mat2 &m) const {
      return {d[0] - m.d[0], d[1] - m.d[1], d[2] - m.d[2], d[3] - m.d[3]};
   }

   Mat2 &operator+=(const Mat2 &m) {
      d[0] += m.d[0];
      d[1] += m.d[1];
      d[2] += m.d[2];
      d[3] += m.d[3];
      return *this;
   }
   Mat2 &operator-=(const Mat2 &m) {
      d[0] -= m.d[0];
      d[1] -= m.d[1];
      d[2] -= m.d[2];
      d[3] -= m.d[3];
      return *this;
   }

   Mat2 operator*(const Mat2 &m) { return Mat2((*this) * m.getCol(0), (*this) * m.getCol(1)); }
   Mat2 &operator*=(const Mat2 &m) {
      this->setCol(0, (*this) * m.getCol(0));
      this->setCol(1, (*this) * m.getCol(1));
      return *this;
   }

   // Comparison
   bool operator==(const Mat2 &m) const {
      return almostEqual(d[0], m.d[0]) && almostEqual(d[1], m.d[1]) && almostEqual(d[2], m.d[2]) &&
             almostEqual(d[3], m.d[3]);
   }
   bool operator!=(const Mat2 &m) const {
      return !almostEqual(d[0], m.d[0]) || !almostEqual(d[1], m.d[1]) ||
             !almostEqual(d[2], m.d[2]) || !almostEqual(d[3], m.d[3]);
   }
};

}  // namespace amath

#endif /* __MAT2_H__ */
