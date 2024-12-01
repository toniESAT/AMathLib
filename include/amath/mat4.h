#ifndef __MAT4_H__
#define __MAT4_H__

#include <math.h>
#include <stdio.h>
#include <vector>

#include "utils.h"

#include "mat3.h"
#include "vec3.h"
#include "vec4.h"

namespace amath {

struct Mat4 {
   scalar d[16];

   /*****************************
    ****  Mat4 constructors  ****
    *****************************/

   Mat4(scalar m0, scalar m1, scalar m2, scalar m3, scalar m4, scalar m5, scalar m6, scalar m7,
        scalar m8, scalar m9, scalar m10, scalar m11, scalar m12, scalar m13, scalar m14,
        scalar m15)
       : d{m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15} {};
   Mat4(scalar k) : Mat4(k, k, k, k, k, k, k, k, k, k, k, k, k, k, k, k) {}
   Mat4() : d{0} {};
   Mat4(const Vec4 &v0, const Vec4 &v1, const Vec4 &v2, const Vec4 &v3) {
      setCol(0, v0);
      setCol(1, v1);
      setCol(2, v2);
      setCol(3, v3);
   }

   static Mat4 identity() { return {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}; }
   static Mat4 nan() { return Mat4(nanf("")); };

   static Mat4 scaling(scalar sx, scalar sy, scalar sz) {
      return {sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0, 0, 0, 0, 1};
   }
   static Mat4 translation(scalar tx, scalar ty, scalar tz) {
      return {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, tx, ty, tz, 1};
   }
   static Mat4 rotationX(scalar a) {
      return {1, 0, 0, 0, 0, cosf(a), -sinf(a), 0, 0, sinf(a), cosf(a), 0, 0, 0, 0, 1};
   }
   static Mat4 rotationY(scalar a) {
      return {cosf(a), 0, sinf(a), 0, 0, 1, 0, 0, -sinf(a), 0, cosf(a), 0, 0, 0, 0, 1};
   }
   static Mat4 rotationZ(scalar a) {
      return {cosf(a), sinf(a), 0, 0, -sinf(a), cosf(a), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
   }

   static Mat4 rotation(scalar x, scalar y, scalar z) {
      // return rotationZ(z) * (rotationY(y) * rotationX(x));

      // Unrolled
      scalar cx = cosf(x);
      scalar sx = sinf(x);
      scalar cy = cosf(y);
      scalar sy = sinf(y);
      scalar cz = cosf(z);
      scalar sz = sinf(z);

      return {
          cx * cy,
          cx * sy * sz - sx * cz,
          cx * sy * cz + sx * sz,
          0,
          sx * cy,
          sx * sy * sz + cx * cz,
          sx * sy * cz - cx * sz,
          0,
          -sy,
          cy * sz,
          cy * cz,
          0,
          0,
          0,
          0,
          1,
      };
   }

   static Mat4 rotationAroundAxis(Vec4 axis, scalar angle) {
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

      return {
          c + x * x * c_c,
          xy * c_c - z * s,
          xz * c_c + y * s,
          0,
          xy * c_c + z * s,
          c + y * y * c_c,
          yz * c_c - x * s,
          0,
          xz * c_c - y * s,
          yz * c_c + x * s,
          c + z * z * c_c,
          0,
          0,
          0,
          0,
          1,
      };
   }

   static Mat4 transform(Vec3 translate = {0, 0, 0}, Vec3 scale = {1, 1, 1}, Vec3 rot = {0, 0, 0}) {
      return Mat4::translation(translate.x(), translate.y(), translate.z()) *
             (Mat4::rotation(rot.x(), rot.y(), rot.z()) *
              Mat4::scaling(scale.x(), scale.y(), scale.z()));
   }

   static Mat4 perspective(scalar fov = PI / 2, scalar aspect = 1, scalar zNear = 1,
                           scalar zFar = 100) {
      return {1.f / (aspect * tanf(fov / 2)),
              0,
              0,
              0,
              0,
              1.f / tanf(fov / 2),
              0,
              0,
              0,
              0,
              zFar / (zFar - zNear),
              1,
              0,
              0,
              (-zNear * zFar) / (zFar - zNear),
              0};
   }

   /*****************************
    *****  Getters/Setters  *****
    *****************************/

   Vec4 getCol(size_t j) const { return {d[4 * j], d[4 * j + 1], d[4 * j + 2], d[4 * j + 3]}; }
   Vec4 getRow(size_t i) const { return {d[i], d[i + 4], d[i + 8], d[i + 12]}; }
   void setCol(size_t j, const Vec4 &v) {
      d[4 * j] = v[0];
      d[4 * j + 1] = v[1];
      d[4 * j + 2] = v[2];
      d[4 * j + 3] = v[3];
   }
   void setRow(size_t i, const Vec4 &v) {
      d[i] = v[0];
      d[i + 4] = v[1];
      d[i + 8] = v[2];
      d[i + 12] = v[3];
   }

   /*****************************
    *********  Methods  *********
    *****************************/

   static size_t size() { return 16; }

   void print() const {
      printf("Mat4\n");
      for (size_t i = 0; i < 4; i++)
         printf(
             "[%.4f][%.4f][%.4f][%.4f]\n", d[4 * 0 + i], d[4 * 1 + i], d[4 * 2 + i], d[4 * 3 + i]);
      printf("\n");
   }

   Mat3 getAdjugate(size_t row, size_t col) const {
      Mat3 adj = Mat3::nan();
      if (col > 3 || row > 3) {
         printf("ERROR at Mat3::getAdjugate: Wrong column or row index");
         return adj;
      };

      size_t adj_idx = 0;
      for (size_t i = 0; i < 16; i++) {
         if (i % 4 != row && i / 4 != col) adj[adj_idx++] = d[i];
      }
      return adj;
   }

   // TODO: unroll this function
   scalar det() const {
      Vec4 sign(1, -1, 1, -1);
      Mat3 adj;
      scalar det = 0;
      // Get det by cofactor expansion for the first row
      for (size_t i = 0; i < 4; i++) {
         adj = getAdjugate(0, i);
         det += sign[i] * d[i * 4] * adj.det();
      }

      return det;
   }

   Mat4 transposed() const {
      return {d[0],
              d[4],
              d[8],
              d[12],
              d[1],
              d[5],
              d[9],
              d[13],
              d[2],
              d[6],
              d[10],
              d[14],
              d[3],
              d[7],
              d[11],
              d[15]};
   };

   std::vector<Vec4> transformPoints(const std::vector<Vec4> &points) const {
      std::vector<Vec4> transformed_points;
      transformed_points.reserve(points.size());
      for (auto p : points) transformed_points.push_back(*this * p);
      return transformed_points;
   }
   std::vector<Vec3> transformPoints(const std::vector<Vec3> &points) const {
      std::vector<Vec3> transformed_points;
      transformed_points.reserve(points.size());
      for (auto p : points) transformed_points.push_back(*this * p);
      return transformed_points;
   }

   /*****************************
    *  Mat4 operator overloads  *
    *****************************/

   // Access
   scalar operator[](size_t i) const { return d[i]; }
   scalar &operator[](size_t i) { return d[i]; }

   scalar operator()(size_t i, size_t j) const {
      if (i >= 0 && i < 4 && j >= 0 && j < 4) return d[i * 4 + j];
   }
   scalar &operator()(size_t i, size_t j) {
      if (i >= 0 && i < 4 && j >= 0 && j < 4) return d[i * 4 + j];
   }

   // Negate
   Mat4 operator-() const {
      return {-d[0],
              -d[1],
              -d[2],
              -d[3],
              -d[4],
              -d[5],
              -d[6],
              -d[7],
              -d[8],
              -d[9],
              -d[10],
              -d[11],
              -d[12],
              -d[13],
              -d[14],
              -d[15]};
   }

   // Operations with scalars
   Mat4 operator+(scalar k) const {
      return {d[0] + k,
              d[1] + k,
              d[2] + k,
              d[3] + k,
              d[4] + k,
              d[5] + k,
              d[6] + k,
              d[7] + k,
              d[8] + k,
              d[9] + k,
              d[10] + k,
              d[11] + k,
              d[12] + k,
              d[13] + k,
              d[14] + k,
              d[15] + k};
   }
   Mat4 operator-(scalar k) const { return *this + (-k); }

   Mat4 &operator+=(scalar k) {
      d[0] += k;
      d[1] += k;
      d[2] += k;
      d[3] += k;
      d[4] += k;
      d[5] += k;
      d[6] += k;
      d[7] += k;
      d[8] += k;
      d[9] += k;
      d[10] += k;
      d[11] += k;
      d[12] += k;
      d[13] += k;
      d[14] += k;
      d[15] += k;
      return *this;
   }
   Mat4 &operator-=(scalar k) {
      *this += (-k);
      return *this;
   }

   Mat4 operator*(scalar k) const {
      return {d[0] * k,
              d[1] * k,
              d[2] * k,
              d[3] * k,
              d[4] * k,
              d[5] * k,
              d[6] * k,
              d[7] * k,
              d[8] * k,
              d[9] * k,
              d[10] * k,
              d[11] * k,
              d[12] * k,
              d[13] * k,
              d[14] * k,
              d[15] * k};
   }
   Mat4 operator/(scalar k) const { return *this * (1 / k); }

   Mat4 operator*=(scalar k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      d[4] *= k;
      d[5] *= k;
      d[6] *= k;
      d[7] *= k;
      d[8] *= k;
      d[9] *= k;
      d[10] *= k;
      d[11] *= k;
      d[12] *= k;
      d[13] *= k;
      d[14] *= k;
      d[15] *= k;
      return *this;
   }
   Mat4 operator/=(scalar k) {
      *this *= (1 / k);
      return *this;
   }

   // Operations with vectors
   Vec4 operator*(const Vec4 &v) const {
      scalar w = d[3] * v.x() + d[7] * v.y() + d[11] * v.z() + d[15];

      return {(d[0] * v.x() + d[4] * v.y() + d[8] * v.z() + d[12]) / w,
              (d[1] * v.x() + d[5] * v.y() + d[9] * v.z() + d[13]) / w,
              (d[2] * v.x() + d[6] * v.y() + d[10] * v.z() + d[14]) / w,
              1};
   }

   Vec3 operator*(const Vec3 &v) const {
      scalar w = d[3] * v.x() + d[7] * v.y() + d[11] * v.z() + d[15];

      return {(d[0] * v.x() + d[4] * v.y() + d[8] * v.z() + d[12]) / w,
              (d[1] * v.x() + d[5] * v.y() + d[9] * v.z() + d[13]) / w,
              (d[2] * v.x() + d[6] * v.y() + d[10] * v.z() + d[14]) / w};
   }

   // Operations with other matrices

   Mat4 operator+(const Mat4 &m) const {
      return {d[0] + m.d[0],
              d[1] + m.d[1],
              d[2] + m.d[2],
              d[3] + m.d[3],
              d[4] + m.d[4],
              d[5] + m.d[5],
              d[6] + m.d[6],
              d[7] + m.d[7],
              d[8] + m.d[8],
              d[9] + m.d[9],
              d[10] + m.d[10],
              d[11] + m.d[11],
              d[12] + m.d[12],
              d[13] + m.d[13],
              d[14] + m.d[14],
              d[15] + m.d[15]};
   }
   Mat4 operator-(const Mat4 &m) const {
      return {d[0] - m.d[0],
              d[1] - m.d[1],
              d[2] - m.d[2],
              d[3] - m.d[3],
              d[4] - m.d[4],
              d[5] - m.d[5],
              d[6] - m.d[6],
              d[7] - m.d[7],
              d[8] - m.d[8],
              d[9] - m.d[9],
              d[10] - m.d[10],
              d[11] - m.d[11],
              d[12] - m.d[12],
              d[13] - m.d[13],
              d[14] - m.d[14],
              d[15] - m.d[15]};
   }

   Mat4 &operator+=(const Mat4 &m) {
      d[0] += m.d[0];
      d[1] += m.d[1];
      d[2] += m.d[2];
      d[3] += m.d[3];
      d[4] += m.d[4];
      d[5] += m.d[5];
      d[6] += m.d[6];
      d[7] += m.d[7];
      d[8] += m.d[8];
      d[9] += m.d[9];
      d[10] += m.d[10];
      d[11] += m.d[11];
      d[12] += m.d[12];
      d[13] += m.d[13];
      d[14] += m.d[14];
      d[15] += m.d[15];
      return *this;
   }
   Mat4 &operator-=(const Mat4 &m) {
      d[0] -= m.d[0];
      d[1] -= m.d[1];
      d[2] -= m.d[2];
      d[3] -= m.d[3];
      d[4] -= m.d[4];
      d[5] -= m.d[5];
      d[6] -= m.d[6];
      d[7] -= m.d[7];
      d[8] -= m.d[8];
      d[9] -= m.d[9];
      d[10] -= m.d[10];
      d[11] -= m.d[11];
      d[12] -= m.d[12];
      d[13] -= m.d[13];
      d[14] -= m.d[14];
      d[15] -= m.d[15];
      return *this;
   }

   Mat4 operator*(const Mat4 &m) {
      Mat4((*this) * m.getCol(0),
           (*this) * m.getCol(1),
           (*this) * m.getCol(2),
           (*this) * m.getCol(3));
   }
   Mat4 &operator*=(const Mat4 &m) {
      this->setCol(0, (*this) * m.getCol(0));
      this->setCol(1, (*this) * m.getCol(1));
      this->setCol(2, (*this) * m.getCol(2));
      this->setCol(3, (*this) * m.getCol(3));
      return *this;
   }

   // Comparison
   // TODO: could this be sped up with SIMD?
   bool operator==(const Mat3 &m) const {
      return almostEqual(d[0], m.d[0]) && almostEqual(d[1], m.d[1]) && almostEqual(d[2], m.d[2]) &&
             almostEqual(d[3], m.d[3]) && almostEqual(d[4], m.d[4]) && almostEqual(d[5], m.d[5]) &&
             almostEqual(d[6], m.d[6]) && almostEqual(d[7], m.d[7]) && almostEqual(d[8], m.d[8]) &&
             almostEqual(d[9], m.d[9]) && almostEqual(d[10], m.d[10]) &&
             almostEqual(d[11], m.d[11]) && almostEqual(d[12], m.d[12]) &&
             almostEqual(d[13], m.d[13]) && almostEqual(d[14], m.d[14]) &&
             almostEqual(d[15], m.d[15]);
   }
   bool operator!=(const Mat3 &m) const {
      return !almostEqual(d[0], m.d[0]) || !almostEqual(d[1], m.d[1]) ||
             !almostEqual(d[2], m.d[2]) || !almostEqual(d[3], m.d[3]) ||
             !almostEqual(d[4], m.d[4]) || !almostEqual(d[5], m.d[5]) ||
             !almostEqual(d[6], m.d[6]) || !almostEqual(d[7], m.d[7]) ||
             !almostEqual(d[8], m.d[8]) || !almostEqual(d[9], m.d[9]) ||
             !almostEqual(d[10], m.d[10]) || !almostEqual(d[11], m.d[11]) ||
             !almostEqual(d[12], m.d[12]) || !almostEqual(d[13], m.d[13]) ||
             !almostEqual(d[14], m.d[14]) || !almostEqual(d[15], m.d[15]);
   }
};

}  // namespace amath

#endif /* __MAT4_H__ */
