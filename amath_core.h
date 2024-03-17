#pragma once

#include <cmath>
#include <fmt/format.h>
#include <immintrin.h>

#define SMALL_NUMBER 1.e-8f

namespace amth {

/********************************/
/****         VECTORS        ****/
/********************************/

struct Vec2 {

   float d[2];

   Vec2() : d{0, 0} {};
   Vec2(float v0, float v1) : d{v0, v1} {};

   float x() const { return d[0]; }
   float y() const { return d[1]; }
   float &x() { return d[0]; }
   float &y() { return d[1]; }

   float squared_length() const { return x() * x() + y() * y(); }
   float length() const { return sqrtf(squared_length()); }

   Vec2 normalized(float tolerance = SMALL_NUMBER) const {
      float l = length();
      if (l > SMALL_NUMBER) return {x() / l, y() / l};
      else return {0, 0};
   }

   bool is_normalized(float tolerance = SMALL_NUMBER) const {
      return squared_length() < SMALL_NUMBER;
   }

   Vec2 perpendicular() const { return {-y(), x()}; }

   void print() const { fmt::print("Vector2 [{},{}]", x(), y()); }

   Vec2 operator-() const { return {-d[0], -d[1]}; }

   Vec2 operator+(Vec2 &v) const { return {x() + v.x(), y() + v.y()}; }
   Vec2 operator-(Vec2 &v) const { return {x() - v.x(), y() - v.y()}; }
   Vec2 operator+(float k) const { return {x() + k, y() + k}; }
   Vec2 operator-(float k) const { return {x() - k, y() - k}; }

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

   Vec2 operator*(float k) const { return {x() * k, y() * k}; }
   Vec2 &operator*=(float k) {
      d[0] *= k;
      d[1] *= k;
      return *this;
   }
   Vec2 operator/(float k) const { return {x() / k, y() / k}; }

   bool operator==(Vec2 &v) const { return (x() == v.x() && y() == v.y()); }
   bool operator!=(Vec2 &v) const { return (x() != v.x() || y() != v.y()); }

   float operator[](int i) const { return d[i]; }
   float &operator[](int i) { return d[i]; }
};

struct Vec3 {

   float d[3];

   Vec3() : d{0, 0, 0} {};
   Vec3(float v0, float v1, float v2) : d{v0, v1, v2} {};

   float x() const { return d[0]; }
   float y() const { return d[1]; }
   float z() const { return d[2]; }
   float &x() { return d[0]; }
   float &y() { return d[1]; }
   float &z() { return d[2]; }

   float squared_length() const { return x() * x() + y() * y() + z() * z(); }
   float length() const { return sqrtf(squared_length()); }

   Vec3 normalized(float tolerance = SMALL_NUMBER) const {
      float l = length();
      if (l > SMALL_NUMBER) return {x() / l, y() / l, z() / l};
      else return {0, 0, 0};
   }

   bool is_normalized(float tolerance = SMALL_NUMBER) const {
      return squared_length() < SMALL_NUMBER;
   }

   void print() const { fmt::print("Vector3 [{},{}, {}]", x(), y(), z()); }

   Vec3 operator-() const { return {-d[0], -d[1], -d[2]}; }

   Vec3 operator-(Vec3 &v) const { return {x() - v.x(), y() - v.y(), z() - v.z()}; }
   Vec3 operator+(Vec3 &v) const { return {x() + v.x(), y() + v.y(), z() + v.z()}; }
   Vec3 operator+(float k) const { return {x() + k, y() + k, z() + k}; }
   Vec3 operator-(float k) const { return {x() - k, y() - k, z() - k}; }

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

   Vec3 operator*(float k) const { return {x() * k, y() * k, z() * k}; }
   Vec3 &operator*=(float k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      return *this;
   }
   Vec3 operator/(float k) const { return {x() / k, y() / k, z() / k}; }

   bool operator==(Vec3 &v) const { return (x() == v.x() && y() == v.y() && z() == v.z()); }
   bool operator!=(Vec3 &v) const { return (x() != v.x() || y() != v.y() || z() != v.z()); }

   float operator[](int i) const { return d[i]; }
   float &operator[](int i) { return d[i]; }
};

struct Vec4 {

   float d[4];

   Vec4() : d{0, 0, 0, 0} {};
   Vec4(float v0, float v1, float v2, float v3) : d{v0, v1, v2, v3} {};

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

   Vec4 normalized(float tolerance = SMALL_NUMBER) const {
      float l = length();
      if (l > SMALL_NUMBER) return {x() / l, y() / l, z() / l, w() / l};
      else return {0, 0, 0, 0};
   }

   bool is_normalized(float tolerance = SMALL_NUMBER) { return squared_length() < SMALL_NUMBER; }

   void print() const { fmt::print("Vector4 [{},{}, {}, {}]", x(), y(), z(), w()); }

   Vec4 operator-() const { return {-d[0], -d[1], -d[2], -d[3]}; }

   Vec4 operator-(Vec4 &v) const { return {x() - v.x(), y() - v.y(), z() - v.z(), w() + v.w()}; }
   Vec4 operator+(Vec4 &v) const { return {x() + v.x(), y() + v.y(), z() + v.z(), w() - v.w()}; }
   Vec4 operator+(float k) const { return {x() + k, y() + k, z() + k, w() + k}; }
   Vec4 operator-(float k) const { return {x() - k, y() - k, z() - k, w() - k}; }

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

   Vec4 operator*(float k) const { return {x() * k, y() * k, z() * k, w() * k}; }
   Vec4 &operator*=(float k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      return *this;
   }
   Vec4 operator/(float k) const { return {x() / k, y() / k, z() / k, w() / k}; }

   bool operator==(Vec4 &v) const {
      return (x() == v.x() && y() == v.y() && z() == v.z() && w() == v.w());
   }
   bool operator!=(Vec4 &v) const {
      return (x() != v.x() || y() != v.y() || z() != v.z() || w() != v.w());
   }

   float operator[](int i) const { return d[i]; }
   float &operator[](int i) { return d[i]; }
};

Vec2 vec_add(const Vec2 &v1, const Vec2 &v2) { return {v1.x() + v2.x(), v1.y() + v2.y()}; }
Vec3 vec_add(const Vec3 &v1, const Vec3 &v2) {
   return {v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z()};
}
Vec4 vec_add(const Vec4 &v1, const Vec4 &v2) {
   return {v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z(), v1.w() + v2.w()};
}

Vec2 vec_mul(const Vec2 &v, float k) { return {v.x() * k, v.y() * k}; }
Vec3 vec_mul(const Vec3 &v, float k) { return {v.x() * k, v.y() * k, v.z() * k}; }
Vec4 vec_mul(const Vec4 &v, float k) { return {v.x() * k, v.y() * k, v.z() * k, v.w() * k}; }

float dot_product(const Vec2 &v1, const Vec2 &v2) { return v1.x() * v2.x() + v1.y() * v2.y(); }
float dot_product(const Vec3 &v1, const Vec3 &v2) {
   return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}
float dot_product(const Vec4 &v1, const Vec4 &v2) {
   return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z() + v1.w() * v2.w();
}

/********************************/
/****        MATRICES        ****/
/********************************/

struct Mat2 {
   float d[4];

   Mat2() : d{0, 0, 0, 0} {};
   Mat2(float m0, float m1, float m2, float m3) : d{m0, m1, m2, m3} {};

   static Mat2 identity() { return {1, 0, 0, 1}; }
   static int size() { return 4; }

   float determinant() const { return d[0] * d[3] - d[1] * d[2]; }

   Vec2 get_column(int j) const { return {d[2 * j], d[2 * j + 1]}; }
   Vec2 get_row(int i) const { return {d[i], d[i + 2]}; }
   void set_column(int j, const Vec2 &v) {
      d[2 * j] = v[0];
      d[2 * j + 1] = v[1];
   }
   void set_row(int i, const Vec2 &v) {
      d[i] = v[0];
      d[i + 2] = v[1];
   }

   Mat2(Vec2 v0, Vec2 v1) {
      set_column(0, v0);
      set_column(1, v1);
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

   Mat2 operator*(float k) const { return {d[0] * k, d[1] * k, d[2] * k, d[3] * k}; }
   Mat2 &operator*=(float k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      return *this;
   }
   Mat2 operator/(float k) const { return *this * (1 / k); }
};

struct Mat3 {
   float d[9];

   Mat3() : d{0, 0, 0, 0, 0, 0, 0, 0, 0} {};
   Mat3(float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8)
       : d{m0, m1, m2, m3, m4, m5, m6, m7, m8} {};

   static Mat3 identity() { return {1, 0, 0, 0, 1, 0, 0, 0, 1}; }

   static Mat3 scaling(float sx, float sy) { return {sx, 0, 0, 0, sy, 0, 0, 0, 1}; }
   static Mat3 translation(float tx, float ty) { return {1, 0, 0, 0, 1, 0, tx, ty, 1}; }
   static Mat3 rotation(float a) { return {cosf(a), -sinf(a), 0, sinf(a), cosf(a), 0, 0, 0, 1}; }

   // TODO
   float determinant() { return d[0] * d[3] - d[1] * d[2]; }

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

   Mat3(Vec3 v0, Vec3 v1, Vec3 v2) {
      set_column(0, v0);
      set_column(1, v1);
      set_column(2, v2);
   }

   float operator[](int i) const { return d[i]; }
   float &operator[](int i) { return d[i]; }

   float operator()(int i, int j) const {
      if (i >= 0 && i < 3 && j >= 0 && j < 3) return d[i * 3 + j];
   }
   float &operator()(int i, int j) {
      if (i >= 0 && i < 3 && j >= 0 && j < 3) return d[i * 3 + j];
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
   Mat3 operator/(float k) const { return *this * (1 / k); }
};

struct Mat4 {
   float d[16];

   Mat4() : d{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} {};
   Mat4(float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7, float m8,
        float m9, float m10, float m11, float m12, float m13, float m14, float m15)
       : d{m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15} {};

   static Mat4 identity() { return {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}; }

   // TODO
   float determinant() const { return d[0] * d[3] - d[1] * d[2]; }

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

   Mat4(Vec4 v0, Vec4 v1, Vec4 v2, Vec4 v3) {
      set_column(0, v0);
      set_column(1, v1);
      set_column(2, v2);
      set_column(3, v3);
   }

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
   Mat4 operator/(float k) const { return *this * (1 / k); }
};

/* Add matrices */
Mat2 mat_add(const Mat2 &m1, const Mat2 &m2) { return m1 + m2; }
Mat3 mat_add(const Mat3 &m1, const Mat3 &m2) { return m1 + m2; }
Mat4 mat_add(const Mat4 &m1, const Mat4 &m2) { return m1 + m2; }

/* Subtract matrices */
Mat2 mat_sub(const Mat2 &m1, const Mat2 &m2) { return m1 - m2; }
Mat3 mat_sub(const Mat3 &m1, const Mat3 &m2) { return m1 - m2; }
Mat4 mat_sub(const Mat4 &m1, const Mat4 &m2) { return m1 - m2; }

/* Multiply matrix by constant */
Mat2 mat_mul(const Mat2 &m1, float k) { return m1 * k; }
Mat3 mat_mul(const Mat3 &m1, float k) { return m1 * k; }
Mat4 mat_mul(const Mat4 &m1, float k) { return m1 * k; }

/* Multiply matrix by vector */
Vec2 mat_mul(const Mat2 &m, const Vec2 &v) {
   // return {m[0] * v[0] + m[2] * v[1], m[1] * v[0] + m[3] * v[1]};
   return {
       dot_product(m.get_row(0), v),
       dot_product(m.get_row(1), v),
   };
}
Vec3 mat_mul(const Mat3 &m, const Vec3 &v) {
   return {
       dot_product(m.get_row(0), v),
       dot_product(m.get_row(1), v),
       dot_product(m.get_row(2), v),
   };
}

Vec4 mat_mul(const Mat4 &m, const Vec4 &v) {
   return {
       dot_product(m.get_row(0), v),
       dot_product(m.get_row(1), v),
       dot_product(m.get_row(2), v),
       dot_product(m.get_row(3), v),
   };
}

/* Multiply matrix by vector SIMD */

// Vec4 mat_mul(const Mat4 &m, const Vec4 &v) {

//    __m256 sum = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//    __m256 v1 = _mm256_set_ps();
//    __m256 v2 = _mm256_set_ps();
//    _mm_fmadd_ps(v1, v2, sum);
//    v1 = _mm256_load_pd();
//    v2 = _mm256_load_pd();
//    _mm_fmadd_ps(v1, v2, sum);
//    v1 = _mm256_load_pd();
//    v2 = _mm256_load_pd();
//    _mm_fmadd_ps(v1, v2, sum);
//    v1 = _mm256_load_pd();
//    v2 = _mm256_load_pd();
//    _mm_fmadd_ps(v1, v2, sum);
// }

/* Multiply matrices */
Mat2 mat_mul(const Mat2 &m1, const Mat2 &m2) {
   Mat2 m3;
   m3.set_column(0, mat_mul(m1, m2.get_column(0)));
   m3.set_column(1, mat_mul(m1, m2.get_column(1)));
   return m3;
}
Mat2 operator*(const Mat2 &m1, const Mat2 &m2) { return mat_mul(m1, m2); }

Mat3 mat_mul(const Mat3 &m1, const Mat3 &m2) {
   Mat3 m3;
   m3.set_column(0, mat_mul(m1, m2.get_column(0)));
   m3.set_column(1, mat_mul(m1, m2.get_column(1)));
   m3.set_column(2, mat_mul(m1, m2.get_column(2)));
   return m3;
}
Mat3 operator*(const Mat3 &m1, const Mat3 &m2) { return mat_mul(m1, m2); }

Mat4 mat_mul(const Mat4 &m1, const Mat4 &m2) {
   Mat4 m3;
   m3.set_column(0, mat_mul(m1, m2.get_column(0)));
   m3.set_column(1, mat_mul(m1, m2.get_column(1)));
   m3.set_column(2, mat_mul(m1, m2.get_column(2)));
   m3.set_column(3, mat_mul(m1, m2.get_column(3)));
   return m3;
}
Mat4 operator*(const Mat4 &m1, const Mat4 &m2) { return mat_mul(m1, m2); }

} // namespace amth