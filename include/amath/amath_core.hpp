#pragma once

#include <math.h>
#include <float.h>
#include <immintrin.h>
#include <vector>

#include "amath_utils.hpp"


namespace amath {

/********************************/
/****         VECTORS        ****/
/********************************/

struct Vec2 {
   float d[2];

   Vec2(float v0, float v1) : d{v0, v1} {};
   Vec2(float v) : d{v, v} {};
   Vec2() : Vec2(0, 0) {};

   float x() const { return d[0]; }
   float y() const { return d[1]; }
   float &x() { return d[0]; }
   float &y() { return d[1]; }

   float squared_length() const { return x() * x() + y() * y(); }
   float length() const { return sqrtf(squared_length()); }

   static Vec2 nan() { return {nanf(""), nanf("")}; }

   Vec2 normalized(float tolerance = FLT_EPSILON) const {
      float l = length();
      if (!is_almost_zero(l, tolerance)) return {x() / l, y() / l};
      else return {0, 0};
   }

   bool is_normalized(const float tolerance = FLT_EPSILON) const {
      return is_almost_zero(squared_length() - 1, tolerance);
   }

   Vec2 perpendicular() const { return {-y(), x()}; }

   void print() const { printf("Vector2 [%.4f, %.4f]\n", x(), y()); }

   /*****************************
    *  Vec2 operator overloads  *
    *****************************/

   // Negate
   Vec2 operator-() const { return {-d[0], -d[1]}; }

   // Operations with scalars
   Vec2 operator+(const float k) const { return {x() + k, y() + k}; }
   Vec2 operator-(const float k) const { return {x() - k, y() - k}; }
   Vec2 operator*(const float k) const { return {x() * k, y() * k}; }
   Vec2 &operator*=(const float k) {
      d[0] *= k;
      d[1] *= k;
      return *this;
   }
   Vec2 operator/(const float k) const { return {x() / k, y() / k}; }
   Vec2 &operator/=(const float k) {
      d[0] /= k;
      d[1] /= k;
      return *this;
   }

   // Operations with vectors
   Vec2 operator+(const Vec2 &v) const { return {x() + v.x(), y() + v.y()}; }
   Vec2 operator-(const Vec2 &v) const { return {x() - v.x(), y() - v.y()}; }
   Vec2 &operator+=(Vec2 &v) {
      x() = x() + v.x();
      y() = y() + v.y();
      return (*this);
   }
   Vec2 &operator-=(Vec2 &v) {
      x() = x() - v.x();
      y() = y() - v.y();
      return (*this);
   }

   // Comparison
   bool operator==(const Vec2 &v) const {
      return fabs(x() - v.x()) < FLT_EPSILON && fabs(y() - v.y()) < FLT_EPSILON;
   }
   bool operator!=(const Vec2 &v) const {
      return fabs(x() - v.x()) > FLT_EPSILON || fabs(y() - v.y()) > FLT_EPSILON;
   }

   // Access
   float operator[](const int i) const { return d[i]; }
   float &operator[](const int i) { return d[i]; }
};

struct Vec3 {
   float d[3];

   Vec3(float v0, float v1, float v2) : d{v0, v1, v2} {};
   Vec3(float v) : d{v, v, v} {};
   Vec3() : Vec3(0, 0, 0) {};

   static Vec3 nan() { return {nanf(""), nanf(""), nanf("")}; }

   static Vec3 up() { return {0, 1, 0}; }

   float x() const { return d[0]; }
   float y() const { return d[1]; }
   float z() const { return d[2]; }
   float &x() { return d[0]; }
   float &y() { return d[1]; }
   float &z() { return d[2]; }

   float squared_length() const { return x() * x() + y() * y() + z() * z(); }
   float length() const { return sqrtf(squared_length()); }

   Vec3 normalized(float tolerance = FLT_EPSILON) const {
      float l = length();
      if (!is_almost_zero(l, tolerance)) return {x() / l, y() / l, z() / l};
      else return {0, 0, 0};
   }

   bool is_normalized(const float tolerance = FLT_EPSILON) const {
      return is_almost_zero(squared_length() - 1, tolerance);
   }

   void print() const { printf("Vector3 [%.4f,%.4f, %.4f]\n", x(), y(), z()); }

   /*****************************
    *  Vec3 operator overloads  *
    *****************************/
   // Negate
   Vec3 operator-() const { return {-d[0], -d[1], -d[2]}; }

   // Operations with scalars
   Vec3 operator+(const float k) const { return {x() + k, y() + k, z() + k}; }
   Vec3 operator-(const float k) const { return {x() - k, y() - k, z() - k}; }
   Vec3 operator*(const float k) const { return {x() * k, y() * k, z() * k}; }
   Vec3 &operator*=(const float k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      return *this;
   }
   Vec3 operator/(const float k) const { return {x() / k, y() / k, z() / k}; }
   Vec3 &operator/=(const float k) {
      d[0] /= k;
      d[1] /= k;
      d[2] /= k;
      return *this;
   }

   // Operations with vectors
   Vec3 operator+(const Vec3 &v) const { return {x() + v.x(), y() + v.y(), z() + v.z()}; }
   Vec3 operator-(const Vec3 &v) const { return {x() - v.x(), y() - v.y(), z() - v.z()}; }
   Vec3 &operator+=(Vec3 &v) {
      x() = x() + v.x();
      y() = y() + v.y();
      z() = z() + v.y();
      return (*this);
   }
   Vec3 &operator-=(Vec3 &v) {
      x() = x() - v.x();
      y() = y() - v.y();
      z() = z() - v.y();
      return (*this);
   }

   // Comparison
   bool operator==(const Vec3 &v) const {
      return fabs(x() - v.x()) < FLT_EPSILON && fabs(y() - v.y()) < FLT_EPSILON &&
             fabs(z() - v.z()) < FLT_EPSILON;
   }
   bool operator!=(const Vec3 &v) const {
      return fabs(x() - v.x()) > FLT_EPSILON || fabs(y() - v.y()) > FLT_EPSILON ||
             fabs(z() - v.z()) > FLT_EPSILON;
   }

   // Access
   float operator[](const int i) const { return d[i]; }
   float &operator[](const int i) { return d[i]; }
};

struct Vec4 {
   float d[4];

   Vec4(float v0, float v1, float v2, float v3) : d{v0, v1, v2, v3} {};
   Vec4(float v) : d{v, v, v, v} {};
   Vec4() : Vec4(0, 0, 0, 0) {};

   static Vec4 nan() { return {nanf(""), nanf(""), nanf(""), nanf("")}; }
   static Vec4 up() { return {0, 1, 0, 0}; }

   float x() const { return d[0]; }
   float y() const { return d[1]; }
   float z() const { return d[2]; }
   float w() const { return d[3]; }
   float &x() { return d[0]; }
   float &y() { return d[1]; }
   float &z() { return d[2]; }
   float &w() { return d[3]; }

   float squared_length() const { return x() * x() + y() * y() + z() * z() + w() * w(); }
   float length() const { return sqrtf(squared_length()); }

   Vec4 normalized(float tolerance = FLT_EPSILON) const {
      float l = length();
      if (!is_almost_zero(l, tolerance)) return {x() / l, y() / l, z() / l, w() / l};
      else return {0, 0, 0, 0};
   }

   Vec4 normalized_homogeneous(float tolerance = FLT_EPSILON) const {
      if (d[3] < tolerance) {
         return normalized();
      } else {
         return *this * (1.f / d[3]);
      }
   }

   bool is_normalized(const float tolerance = FLT_EPSILON) {
      return is_almost_zero(squared_length() - 1, tolerance);
   }

   void print() const { printf("Vector4 [%.4f, %.4f, %.4f, %.4f]\n", x(), y(), z(), w()); }

   float dot_product(const Vec4 &v2) {
      return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z() + v1.w() * v2.w();
   }

   /*****************************
    *  Vec4 operator overloads  *
    *****************************/
   // Negate
   Vec4 operator-() const { return {-d[0], -d[1], -d[2], -d[3]}; }

   // Operations with scalars
   Vec4 operator+(const float k) const { return {x() + k, y() + k, z() + k, w() + k}; }
   Vec4 operator-(const float k) const { return {x() - k, y() - k, z() - k, w() - k}; }
   Vec4 operator*(const float k) const { return {x() * k, y() * k, z() * k, w() * k}; }
   Vec4 &operator*=(const float k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      return *this;
   }
   Vec4 operator/(const float k) const { return {x() / k, y() / k, z() / k, w() / k}; }
   Vec4 &operator/=(const float k) {
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
      return (*this);
   }
   Vec4 &operator-=(Vec4 &v) {
      x() -= v.x();
      y() -= v.y();
      z() -= v.z();
      w() -= v.w();
      return (*this);
   }

   // Comparison
   bool operator==(const Vec4 &v) const {
      return fabs(x() - v.x()) < FLT_EPSILON && fabs(y() - v.y()) < FLT_EPSILON &&
             fabs(z() - v.z()) < FLT_EPSILON && fabs(w() - v.w()) < FLT_EPSILON;
   }
   bool operator!=(const Vec4 &v) const {
      return fabs(x() - v.x()) > FLT_EPSILON || fabs(y() - v.y()) > FLT_EPSILON ||
             fabs(z() - v.z()) > FLT_EPSILON || fabs(w() - v.w()) > FLT_EPSILON;
   }

   // Access
   float operator[](const int i) const { return d[i]; }
   float &operator[](const int i) { return d[i]; }
};

inline float dot_product(const Vec2 &v1, const Vec2 &v2) {
   return v1.x() * v2.x() + v1.y() * v2.y();
}
inline float dot_product(const Vec3 &v1, const Vec3 &v2) {
   return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

/********************************/
/****        MATRICES        ****/
/********************************/

struct Mat2 {
   float d[4];

   Mat2(float m0, float m1, float m2, float m3) : d{m0, m1, m2, m3} {};
   Mat2(float k) : Mat2(k, k, k, k) {}
   Mat2() : d{0} {};
   Mat2(const Vec2 &v0, const Vec2 &v1) {
      set_column(0, v0);
      set_column(1, v1);
   }

   static int size() { return 4; }
   static Mat2 identity() { return {1, 0, 0, 1}; }
   static Mat2 nan() { return Mat2(nanf("")); }

   void print() const { printf("2D Matrix\n[%.4f][%.4f]\n[%.4f][%.4f]", d[0], d[2], d[1], d[3]); }

   float determinant() const { return d[0] * d[3] - d[1] * d[2]; }  // By Sarrus' rule

   Vec2 get_column(const int j) const { return {d[2 * j], d[2 * j + 1]}; }
   Vec2 get_row(const int i) const { return {d[i], d[i + 2]}; }
   void set_column(const int j, const Vec2 &v) {
      d[2 * j] = v[0];
      d[2 * j + 1] = v[1];
   }
   void set_row(const int i, const Vec2 &v) {
      d[i] = v[0];
      d[i + 2] = v[1];
   }

   float operator[](int i) const { return d[i]; }
   float &operator[](int i) { return d[i]; }

   float operator()(int i, int j) const {
      if (i >= 0 && i < 2 && j >= 0 && j < 2) return d[i * 2 + j];
   }
   float &operator()(int i, int j) {
      if (i >= 0 && i < 2 && j >= 0 && j < 2) return d[i * 2 + j];
   }

   Mat2 operator+(const Mat2 &m) const {
      return {d[0] + m.d[0], d[1] + m.d[1], d[2] + m.d[2], d[3] + m.d[3]};
   }
   Mat2 operator-(const Mat2 &m) const {
      return {d[0] - m.d[0], d[1] - m.d[1], d[2] - m.d[2], d[3] - m.d[3]};
   }

   Mat2 operator*(const float k) const { return {d[0] * k, d[1] * k, d[2] * k, d[3] * k}; }
   Mat2 operator*(const Mat2 &m);
   Mat2 &operator*=(const float k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      return *this;
   }
   Mat2 operator/(const float k) const { return *this * (1 / k); }
};

struct Mat3 {
   float d[9];

   Mat3(float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8)
       : d{m0, m1, m2, m3, m4, m5, m6, m7, m8} {};
   Mat3(float k) : Mat3(k, k, k, k, k, k, k, k, k) {}
   Mat3() : d{0} {};
   Mat3(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2) {
      set_column(0, v0);
      set_column(1, v1);
      set_column(2, v2);
   }

   static int size() { return 9; }
   static Mat3 identity() { return {1, 0, 0, 0, 1, 0, 0, 0, 1}; }
   static Mat3 nan() { return Mat3(nanf("")); };

   // 2D transformation matrices (for 2D homogeneous coordinate systems)
   static Mat3 scaling(float sx, float sy) { return {sx, 0, 0, 0, sy, 0, 0, 0, 1}; }
   static Mat3 translation(float tx, float ty) { return {1, 0, 0, 0, 1, 0, tx, ty, 1}; }
   static Mat3 rotation(float a) { return {cosf(a), -sinf(a), 0, sinf(a), cosf(a), 0, 0, 0, 1}; }

   // 3D transformation matrices (for 3D euclidean coordinate systems)
   static Mat3 rotation_around_axis(Vec3 axis, float angle) {
      if (!axis.is_normalized()) axis = axis.normalized();

      // Sine, cosine and complements of angle
      float c = cosf(angle);
      float s = sinf(angle);
      float c_c = 1 - c;
      float s_c = 1 - s;

      // Products
      float x = axis.x();
      float y = axis.y();
      float z = axis.z();
      float xy = x * y;
      float yz = y * z;
      float xz = x * z;

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

   void print() const {
      printf("3D Matrix\n");
      for (int i = 0; i < 3; i++)
         printf("[%.4f][%.4f][%.4f]\n", d[3 * i], d[3 * i + 1], d[3 * i + 2]);
   }

   float determinant() {  // By Sarrus' rule
      return d[0] * d[4] * d[8] + d[2] * d[3] * d[7] + d[1] * d[5] * d[6] - d[2] * d[4] * d[6] -
             d[1] * d[3] * d[8] - d[0] * d[5] * d[7];
   }

   Vec3 get_column(int j) const { return {d[3 * j], d[3 * j + 1], d[3 * j + 2]}; }
   Vec3 get_row(int i) const { return {d[i], d[i + 3], d[i + 6]}; }
   void set_column(int j, const Vec3 &v) {
      d[3 * j] = v[0];
      d[3 * j + 1] = v[1];
      d[3 * j + 2] = v[2];
   }
   void set_row(int i, const Vec3 &v) {
      d[i] = v[0];
      d[i + 3] = v[1];
      d[i + 6] = v[2];
   }

   float operator[](int i) const { return d[i]; }
   float &operator[](int i) { return d[i]; }

   float operator()(size_t i, size_t j) const {
      if (i < 3 && j < 3) return d[i * 3 + j];
      else {
         printf("ERROR at Mat3(): Wrong row or column index.\n");
         return nanf("");
      }
   }
   float &operator()(size_t i, size_t j) {
      if (i > 3 || j < 3) return d[i * 3 + j];
      else {
         printf("ERROR at Mat3(): Wrong row or column index.\n");
         exit(1);  // ! should this throw?
      }
   }

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
   Mat3 operator-(const Mat3 &m) const {
      return {d[0] - m.d[0],
              d[1] - m.d[1],
              d[2] - m.d[2],
              d[3] - m.d[3],
              d[4] - m.d[4],
              d[5] - m.d[5],
              d[6] - m.d[6],
              d[7] - m.d[7],
              d[8] - m.d[8]};
   }

   Mat3 operator*(float k) const {
      return {
          d[0] * k, d[1] * k, d[2] * k, d[3] * k, d[4] * k, d[5] * k, d[6] * k, d[7] * k, d[8] * k};
   }

   Mat3 operator*(const Mat3 &m);
   Vec3 operator*(const Vec3 &v);

   Mat3 &operator*=(float k) {
      for (int i = 0; i < this->size(); i++) d[i] *= k;
      return *this;
   }
   Mat3 operator/(float k) const { return *this * (1 / k); }
};

struct Mat4 {
   float d[16];

   Mat4(float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8,
        float m9, float m10, float m11, float m12, float m13, float m14, float m15)
       : d{m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15} {};
   Mat4(float k) : Mat4(k, k, k, k, k, k, k, k, k, k, k, k, k, k, k, k) {}
   Mat4() : d{0} {};
   Mat4(const Vec4 &v0, const Vec4 &v1, const Vec4 &v2, const Vec4 &v3) {
      set_column(0, v0);
      set_column(1, v1);
      set_column(2, v2);
      set_column(3, v3);
   }

   static int size() { return 16; }
   static Mat4 identity() { return {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}; }
   static Mat4 nan() { return Mat4(nanf("")); };

   static Mat4 scaling(float sx, float sy, float sz) {
      return {sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0, 0, 0, 0, 1};
   }
   static Mat4 translation(float tx, float ty, float tz) {
      return {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, tx, ty, tz, 1};
   }
   static Mat4 rotationX(float a) {
      return {1, 0, 0, 0, 0, cosf(a), -sinf(a), 0, 0, sinf(a), cosf(a), 0, 0, 0, 0, 1};
   }
   static Mat4 rotationY(float a) {
      return {cosf(a), 0, sinf(a), 0, 0, 1, 0, 0, -sinf(a), 0, cosf(a), 0, 0, 0, 0, 1};
   }
   static Mat4 rotationZ(float a) {
      return {cosf(a), sinf(a), 0, 0, -sinf(a), cosf(a), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
   }

   static Mat4 rotation(float x, float y, float z);

   static Mat4 rotation_around_axis(Vec4 axis, float angle) {
      if (!axis.is_normalized()) axis = axis.normalized();

      // Sine, cosine and complements of angle
      float c = cosf(angle);
      float s = sinf(angle);
      float c_c = 1 - c;
      float s_c = 1 - s;

      // Products
      float x = axis.x();
      float y = axis.y();
      float z = axis.z();
      float xy = x * y;
      float yz = y * z;
      float xz = x * z;

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

   static Mat4 transform(Vec3 translate = {0, 0, 0}, Vec3 scale = {1, 1, 1}, Vec3 rot = {0, 0, 0});

   static Mat4 perspective(float fov = PI / 2, float aspect = 1, float zNear = 1,
                           float zFar = 100) {
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

   void print() const {
      printf("4D Matrix\n");
      for (int i = 0; i < 4; i++)
         printf(
             "[%.4f][%.4f][%.4f][%.4f]\n", d[4 * 0 + i], d[4 * 1 + i], d[4 * 2 + i], d[4 * 3 + i]);
      printf("\n");
   }

   Mat3 get_adjugate(size_t row, size_t col) const {
      Mat3 adj = Mat3::nan();
      if (col > 3 || row > 3) {
         printf("ERROR at Mat3::get_adjugate: Wrong column or row index");
         return adj;
      };

      size_t adj_idx = 0;
      for (int i = 0; i < 16; i++) {
         if (i % 4 != row && i / 4 != col) adj[adj_idx++] = d[i];
      }
      return adj;
   }

   float determinant() const {
      Vec4 sign(1, -1, 1, -1);
      Mat3 adj;
      float det = 0;
      // Get det by cofactor expansion for the first row
      for (int i = 0; i < 4; i++) {
         adj = get_adjugate(0, i);
         det += sign[i] * d[i * 4] * adj.determinant();
      }

      return det;
   }

   Vec4 get_column(int j) const { return {d[4 * j], d[4 * j + 1], d[4 * j + 2], d[4 * j + 3]}; }
   Vec4 get_row(int i) const { return {d[i], d[i + 4], d[i + 8], d[i + 12]}; }
   void set_column(int j, const Vec4 &v) {
      d[4 * j] = v[0];
      d[4 * j + 1] = v[1];
      d[4 * j + 2] = v[2];
      d[4 * j + 3] = v[3];
   }
   void set_row(int i, const Vec4 &v) {
      d[i] = v[0];
      d[i + 4] = v[1];
      d[i + 8] = v[2];
      d[i + 12] = v[3];
   }

   std::vector<Vec4> transform_points(const std::vector<Vec4> &points) const;

   // Operator overloads
   float operator[](int i) const { return d[i]; }
   float &operator[](int i) { return d[i]; }

   float operator()(int i, int j) const {
      if (i >= 0 && i < 4 && j >= 0 && j < 4) return d[i * 4 + j];
   }
   float &operator()(int i, int j) {
      if (i >= 0 && i < 4 && j >= 0 && j < 4) return d[i * 4 + j];
   }

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

   Mat4 operator*(float k) const {
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
   Mat4 operator*(const Mat4 &m);
   Vec4 operator*(const Vec4 &v);
   Mat4 &operator*=(float k) {
      for (int i = 0; i < this->size(); i++) d[i] *= k;
      return *this;
   }
   Mat4 operator/(float k) const { return *this * (1 / k); }
};

/* Add matrices */
inline Mat2 mat_add(const Mat2 &m1, const Mat2 &m2) { return m1 + m2; }
inline Mat3 mat_add(const Mat3 &m1, const Mat3 &m2) { return m1 + m2; }
inline Mat4 mat_add(const Mat4 &m1, const Mat4 &m2) { return m1 + m2; }

/* Subtract matrices */
inline Mat2 mat_sub(const Mat2 &m1, const Mat2 &m2) { return m1 - m2; }
inline Mat3 mat_sub(const Mat3 &m1, const Mat3 &m2) { return m1 - m2; }
inline Mat4 mat_sub(const Mat4 &m1, const Mat4 &m2) { return m1 - m2; }

/* Multiply matrix by constant */
inline Mat2 mat_mul(const Mat2 &m1, float k) { return m1 * k; }
inline Mat3 mat_mul(const Mat3 &m1, float k) { return m1 * k; }
inline Mat4 mat_mul(const Mat4 &m1, float k) { return m1 * k; }

/* Multiply matrix by vector */
inline Vec2 mat_mul(const Mat2 &m, const Vec2 &v) {
   // return {m[0] * v[0] + m[2] * v[1], m[1] * v[0] + m[3] * v[1]};
   return {
       dot_product(m.get_row(0), v),
       dot_product(m.get_row(1), v),
   };
}
inline Vec3 mat_mul(const Mat3 &m, const Vec3 &v) {
   return {
       dot_product(m.get_row(0), v),
       dot_product(m.get_row(1), v),
       dot_product(m.get_row(2), v),
   };
}

inline Vec3 Mat3::operator*(const Vec3 &v) { return mat_mul(*this, v); }

inline Vec4 mat_mul(const Mat4 &m, const Vec4 &v) {
   return {
       dot_product(m.get_row(0), v),
       dot_product(m.get_row(1), v),
       dot_product(m.get_row(2), v),
       dot_product(m.get_row(3), v),
   };
}
inline Vec4 Mat4::operator*(const Vec4 &v) { return mat_mul(*this, v); }

/* Multiply matrices */
inline Mat2 mat_mul(const Mat2 &m1, const Mat2 &m2) {
   Mat2 m3;
   m3.set_column(0, mat_mul(m1, m2.get_column(0)));
   m3.set_column(1, mat_mul(m1, m2.get_column(1)));
   return m3;
}
inline Mat2 Mat2::operator*(const Mat2 &m) { return mat_mul(*this, m); }

inline Mat3 mat_mul(const Mat3 &m1, const Mat3 &m2) {
   Mat3 m3;
   m3.set_column(0, mat_mul(m1, m2.get_column(0)));
   m3.set_column(1, mat_mul(m1, m2.get_column(1)));
   m3.set_column(2, mat_mul(m1, m2.get_column(2)));
   return m3;
}

inline Mat3 Mat3::operator*(const Mat3 &m) { return mat_mul(*this, m); }

inline Mat4 mat_mul(const Mat4 &m1, const Mat4 &m2) {
   Mat4 m3;
   m3.set_column(0, mat_mul(m1, m2.get_column(0)));
   m3.set_column(1, mat_mul(m1, m2.get_column(1)));
   m3.set_column(2, mat_mul(m1, m2.get_column(2)));
   m3.set_column(3, mat_mul(m1, m2.get_column(3)));
   return m3;
}

inline Mat4 Mat4::operator*(const Mat4 &m) { return mat_mul(*this, m); }

inline Mat4 mat_concat(std::vector<Mat4> matrices) {
   Mat4 r = Mat4::identity();
   for (size_t i = 0; i < matrices.size(); i++) r = matrices[i] * r;
   return r;
}

// Cross product
inline Vec3 cross_product(const Vec3 &v1, const Vec3 &v2) {
   return {Vec3(v1.y() * v2.z() - v1.z() * v2.y(),
                v1.z() * v2.x() - v1.x() * v2.z(),
                v1.x() * v2.y() - v1.y() * v2.x())};
}

inline Vec4 cross_product(const Vec4 &v1, const Vec4 &v2) {
   return {Vec4(v1.y() * v2.z() - v1.z() * v2.y(),
                v1.z() * v2.x() - v1.x() * v2.z(),
                v1.x() * v2.y() - v1.y() * v2.x(),
                0)};
}

inline std::vector<Vec4> Mat4::transform_points(const std::vector<Vec4> &points) const {
   std::vector<Vec4> transformed_points;
   transformed_points.reserve(points.size());
   for (auto p : points) transformed_points.push_back(mat_mul(*this, p));
   return transformed_points;
}

inline Mat4 Mat4::rotation(float x, float y, float z) {
   return mat_mul(rotationZ(z), mat_mul(rotationY(y), rotationX(x)));

   float cx = cosf(x);
   float sx = sinf(x);
   float cy = cosf(y);
   float sy = sinf(y);
   float cz = cosf(z);
   float sz = sinf(z);

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

inline Mat4 Mat4::transform(Vec3 translate, Vec3 scale, Vec3 rot) {
   Mat4 tr = Mat4::scaling(scale.x(), scale.y(), scale.z());
   tr = mat_mul(Mat4::rotation(rot.x(), rot.y(), rot.z()), tr);
   tr = mat_mul(Mat4::translation(translate.x(), translate.y(), translate.z()), tr);
   return tr;
}

}  // namespace amath