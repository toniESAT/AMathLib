#ifndef __MAT3_H__
#define __MAT3_H__

#include <math.h>
#include <stdio.h>

#include "common_defs.h"

#include "vec3.h"

namespace amath {

struct Mat3 {
   scalar d[9];

   /*****************************
    ****  Mat3 constructors  ****
    *****************************/

   Mat3(scalar m0, scalar m1, scalar m2, scalar m3, scalar m4, scalar m5, scalar m6, scalar m7,
        scalar m8)
       : d{m0, m1, m2, m3, m4, m5, m6, m7, m8} {};
   Mat3(scalar k) : Mat3(k, k, k, k, k, k, k, k, k) {}
   Mat3() : d{0} {};
   Mat3(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2) {
      setCol(0, v0);
      setCol(1, v1);
      setCol(2, v2);
   }

   // 2D transformation matrices (for 2D homogeneous coordinate systems)
   static Mat3 scaling(scalar sx, scalar sy) { return {sx, 0, 0, 0, sy, 0, 0, 0, 1}; }
   static Mat3 translation(scalar tx, scalar ty) { return {1, 0, 0, 0, 1, 0, tx, ty, 1}; }
   static Mat3 rotation(scalar a) { return {cosf(a), -sinf(a), 0, sinf(a), cosf(a), 0, 0, 0, 1}; }

   // 3D transformation matrices (for 3D euclidean coordinate systems)
   static Mat3 rotationAroundAxis(Vec3 axis, scalar angle) {
      if (!axis.isNormalized()) axis = axis.normalized();

      // Sine, cosine and complements of angle
      scalar c = cosf(angle);
      scalar s = sinf(angle);
      scalar c_c = 1 - c;
      scalar s_c = 1 - s;

      // Products
      scalar x = axis.x();
      scalar y = axis.y();
      scalar z = axis.z();
      scalar xy = x * y;
      scalar yz = y * z;
      scalar xz = x * z;

      return {c + x * x * c_c,
              xy * c_c - z * s,
              xz * c_c + y * s,
              xy * c_c + z * s,
              c + y * y * c_c,
              yz * c_c - x * s,
              xz * c_c - y * s,
              yz * c_c + x * s,
              c + z * z * c_c};
   }

   /*****************************
    *****  Getters/Setters  *****
    *****************************/

   Vec3 getCol(int j) const { return {d[3 * j], d[3 * j + 1], d[3 * j + 2]}; }
   Vec3 getRow(int i) const { return {d[i], d[i + 3], d[i + 6]}; }
   void setCol(int j, const Vec3 &v) {
      d[3 * j] = v[0];
      d[3 * j + 1] = v[1];
      d[3 * j + 2] = v[2];
   }
   void setRow(int i, const Vec3 &v) {
      d[i] = v[0];
      d[i + 3] = v[1];
      d[i + 6] = v[2];
   }

   /*****************************
    *********  Methods  *********
    *****************************/

   static int size() { return 9; }
   static Mat3 identity() { return {1, 0, 0, 0, 1, 0, 0, 0, 1}; }
   static Mat3 nan() { return Mat3(nanf("")); };

   scalar det() {  // By Sarrus' rule
      return d[0] * d[4] * d[8] + d[2] * d[3] * d[7] + d[1] * d[5] * d[6] - d[2] * d[4] * d[6] -
             d[1] * d[3] * d[8] - d[0] * d[5] * d[7];
   }

   Mat3 transposed() const { return {d[0], d[3], d[6], d[1], d[4], d[7], d[2], d[5], d[8]}; };

   void print() const {
      printf("Mat3:\n");
      for (int i = 0; i < 3; i++)
         printf("[%.4f][%.4f][%.4f]\n", d[3 * i], d[3 * i + 1], d[3 * i + 2]);
   }

   /*****************************
    *  Mat3 operator overloads  *
    *****************************/

   // Access
   scalar operator[](int i) const { return d[i]; }
   scalar &operator[](int i) { return d[i]; }

   scalar operator()(size_t i, size_t j) const {
      // if (i < 3 && j < 3)
      return d[i * 3 + j];
      // else {
      //    printf("ERROR at Mat3(): Wrong row or column index.\n");
      //    return nanf("");
      // }
   }
   scalar &operator()(size_t i, size_t j) {
      // if (i > 3 || j < 3)
      return d[i * 3 + j];
      // else {
      //    printf("ERROR at Mat3(): Wrong row or column index.\n");
      //    exit(1);
      // }
   }

   // Negate
   Mat3 operator-() const {
      return {-d[0], -d[1], -d[2], -d[3], -d[4], -d[5], -d[6], -d[7], -d[8]};
   }

   // Operations with scalars
   Mat3 operator+(const scalar k) const {
      return {
          d[0] + k, d[1] + k, d[2] + k, d[3] + k, d[4] + k, d[5] + k, d[6] + k, d[7] + k, d[8] + k};
   }
   Mat3 operator-(const scalar k) const { return *this + (-k); }

   Mat3 &operator+=(const scalar k) {
      d[0] += k;
      d[1] += k;
      d[2] += k;
      d[3] += k;
      d[4] += k;
      d[5] += k;
      d[6] += k;
      d[7] += k;
      d[8] += k;
      return *this;
   }
   Mat3 &operator-=(const scalar k) {
      *this += (-k);
      return *this;
   }

   Mat3 operator*(const scalar k) const {
      return {
          d[0] * k, d[1] * k, d[2] * k, d[3] * k, d[4] * k, d[5] * k, d[6] * k, d[7] * k, d[8] * k};
   }
   Mat3 operator/(const scalar k) const { return *this * (1 / k); }
   Mat3 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      d[4] *= k;
      d[5] *= k;
      d[6] *= k;
      d[7] *= k;
      d[8] *= k;
      return *this;
   }
   Mat3 &operator/=(const scalar k) {
      *this *= (1 / k);
      return *this;
   }

   // Operations with vectors
   Vec3 operator*(const Vec3 &v) const {
      return {
          this->getRow(0).dot(v),
          this->getRow(1).dot(v),
          this->getRow(2).dot(v),
      };
   }

   // Operations with other matrices
   Mat3 operator+(const Mat3 &m) const {
      return {d[0] + m.d[0],
              d[1] + m.d[1],
              d[2] + m.d[2],
              d[3] + m.d[3],
              d[4] + m.d[4],
              d[5] + m.d[5],
              d[6] + m.d[6],
              d[7] + m.d[7],
              d[8] + m.d[8]};
   }
   Mat3 operator-(const Mat3 &m) const { return *this + (-m); }

   Mat3 &operator+=(const Mat3 &m) {
      d[0] += m.d[0];
      d[1] += m.d[1];
      d[2] += m.d[2];
      d[3] += m.d[3];
      d[4] += m.d[4];
      d[5] += m.d[5];
      d[6] += m.d[6];
      d[7] += m.d[7];
      d[8] += m.d[8];
      return *this;
   }
   Mat3 &operator-=(const Mat3 &m) {
      *this += (-m);
      return *this;
   }

   Mat3 operator*(const Mat3 &m) const {
      return Mat3((*this) * m.getCol(0), (*this) * m.getCol(1), (*this) * m.getCol(2));
   }
   Mat3 &operator*=(const Mat3 &m) {
      this->setCol(0, (*this) * m.getCol(0));
      this->setCol(1, (*this) * m.getCol(1));
      this->setCol(2, (*this) * m.getCol(2));
      return *this;
   }
};

}  // namespace amath

#endif /* __MAT3_H__ */
