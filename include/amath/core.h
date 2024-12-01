#pragma once

#include <math.h>
#include <vector>

#include "utils.h"

namespace amath {

/********************************/
/****         VECTORS        ****/
/********************************/

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
   scalar operator[](const int i) const { return d[i]; }
   scalar &operator[](const int i) { return d[i]; }
};

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
   scalar operator[](const int i) const { return d[i]; }
   scalar &operator[](const int i) { return d[i]; }
};

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
      if (isAlmostZero(d[3]))
         //  if (isAlmostZero(l, tolerance)) return {0, 0, 0, 0};  // TODO: Prevent div by 0 ???
         return {x() * k, y() * k, z() * k, 0};
      // If homogeneous point (w!=1), normalize to w = 1
      return *this * (1.f / d[3]);
   }

   bool isNormalized(const scalar tolerance = AM_EPSILON) {
      return isAlmostZero(lengthSquared() - 1, tolerance);
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
   scalar operator[](const int i) const { return d[i]; }
   scalar &operator[](const int i) { return d[i]; }
};

/********************************/
/****        MATRICES        ****/
/********************************/

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

   Vec2 getCol(const int j) const { return {d[2 * j], d[2 * j + 1]}; }
   Vec2 getRow(const int i) const { return {d[i], d[i + 2]}; }
   void setCol(const int j, const Vec2 &v) {
      d[2 * j] = v[0];
      d[2 * j + 1] = v[1];
   }
   void setRow(const int i, const Vec2 &v) {
      d[i] = v[0];
      d[i + 2] = v[1];
   }

   /*****************************
    *********  Methods  *********
    *****************************/

   static int size() { return 4; }

   scalar det() const { return d[0] * d[3] - d[1] * d[2]; }

   Mat2 transposed() const { return {d[0], d[2], d[1], d[3]}; };

   void print() const { printf("Mat2:\n[%.4f][%.4f]\n[%.4f][%.4f]", d[0], d[2], d[1], d[3]); }

   /*****************************
    *  Mat2 operator overloads  *
    *****************************/

   // Access
   scalar operator[](int i) const { return d[i]; }
   scalar &operator[](int i) { return d[i]; }

   scalar operator()(int i, int j) const {
      if (i >= 0 && i < 2 && j >= 0 && j < 2) return d[i * 2 + j];
   }
   scalar &operator()(int i, int j) {
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
      return {
          this->getRow(0).dot(v),
          this->getRow(1).dot(v),
      };
   }

   // Operations with other matrices
   Mat2 operator+(const Mat2 &m) const {
      return {d[0] + m.d[0], d[1] + m.d[1], d[2] + m.d[2], d[3] + m.d[3]};
   }
   Mat2 operator-(const Mat2 &m) const { return *this + (-m); }

   Mat2 &operator+=(const Mat2 &m) {
      d[0] += m.d[0];
      d[1] += m.d[1];
      d[2] += m.d[2];
      d[3] += m.d[3];
      return *this;
   }
   Mat2 &operator-=(const Mat2 &m) {
      *this += (-m);
      return *this;
   }

   Mat2 operator*(const Mat2 &m) { return Mat2((*this) * m.getCol(0), (*this) * m.getCol(1)); }
   Mat2 &operator*=(const Mat2 &m) {
      this->setCol(0, (*this) * m.getCol(0));
      this->setCol(1, (*this) * m.getCol(1));
      return *this;
   }
};

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

   Vec4 getCol(int j) const { return {d[4 * j], d[4 * j + 1], d[4 * j + 2], d[4 * j + 3]}; }
   Vec4 getRow(int i) const { return {d[i], d[i + 4], d[i + 8], d[i + 12]}; }
   void setCol(int j, const Vec4 &v) {
      d[4 * j] = v[0];
      d[4 * j + 1] = v[1];
      d[4 * j + 2] = v[2];
      d[4 * j + 3] = v[3];
   }
   void setRow(int i, const Vec4 &v) {
      d[i] = v[0];
      d[i + 4] = v[1];
      d[i + 8] = v[2];
      d[i + 12] = v[3];
   }

   /*****************************
    *********  Methods  *********
    *****************************/

   static int size() { return 16; }

   void print() const {
      printf("Mat4\n");
      for (int i = 0; i < 4; i++)
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
      for (int i = 0; i < 16; i++) {
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
      for (int i = 0; i < 4; i++) {
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
   scalar operator[](int i) const { return d[i]; }
   scalar &operator[](int i) { return d[i]; }

   scalar operator()(int i, int j) const {
      if (i >= 0 && i < 4 && j >= 0 && j < 4) return d[i * 4 + j];
   }
   scalar &operator()(int i, int j) {
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
   // TODO: unroll these functions
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
   Mat4 operator-(const Mat4 &m) const { return *this + (-m); }

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
   Mat4 &operator+=(const Mat4 &m) {
      *this += (-m);
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
};

}  // namespace amath