#ifndef __MAT3_H__
#define __MAT3_H__

#include <math.h>
#include <stdio.h>
#include <vector>

#include "utils.h"

#include "vec3.h"
#include "vec2.h"

namespace amath {

/**
 * @brief A 3x3 matrix class for linear algebra operations and transformations
 *
 * This class represents a 3x3 matrix and provides various operations including:
 * - Basic arithmetic operations (+, -, *, /)
 * - Matrix-vector multiplication
 * - Matrix-matrix multiplication
 * - Determinant calculation
 * - Transpose operation
 * - 2D and 3D point transformations
 * - Special transformation matrices (scaling, rotation, translation)
 */
struct Mat3 {
   scalar d[9];  ///< Array storing matrix elements in column-major order

   /**
    * @name Constructors
    * Different ways to construct a Mat3 object
    */
   ///@{

   /**
    * @brief Constructs a matrix with specified elements
    * @param m0 Element at position (0,0)
    * @param m1 Element at position (1,0)
    * @param m2 Element at position (2,0)
    * @param m3 Element at position (0,1)
    * @param m4 Element at position (1,1)
    * @param m5 Element at position (2,1)
    * @param m6 Element at position (0,2)
    * @param m7 Element at position (1,2)
    * @param m8 Element at position (2,2)
    */
   Mat3(scalar m0, scalar m1, scalar m2, scalar m3, scalar m4, scalar m5, scalar m6, scalar m7,
        scalar m8)
       : d{m0, m1, m2, m3, m4, m5, m6, m7, m8} {};

   /**
    * @brief Constructs a matrix with all elements set to the same value
    * @param k Value to set all elements to
    */
   Mat3(scalar k) : Mat3(k, k, k, k, k, k, k, k, k) {}

   /**
    * @brief Default constructor - initializes all elements to zero
    */
   Mat3() : d{0} {};

   /**
    * @brief Constructs a matrix from three column vectors
    * @param v0 First column vector
    * @param v1 Second column vector
    * @param v2 Third column vector
    */
   Mat3(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2) {
      setCol(0, v0);
      setCol(1, v1);
      setCol(2, v2);
   }
   ///@}

   /**
    * @name Static Factory Methods
    * Methods for creating special matrices
    */
   ///@{

   /**
    * @brief Creates a 3x3 identity matrix
    * @return Identity matrix
    */
   static Mat3 identity() { return {1, 0, 0, 0, 1, 0, 0, 0, 1}; }

   /**
    * @brief Creates a matrix with all elements set to NaN
    * @return Matrix with NaN elements
    */
   static Mat3 nan() { return Mat3(nanf("")); };

   /**
    * @brief Creates a 2D scaling transformation matrix
    * @param sx Scale factor in x direction
    * @param sy Scale factor in y direction
    * @return Scaling transformation matrix
    */
   static Mat3 scaling(scalar sx, scalar sy) { return {sx, 0, 0, 0, sy, 0, 0, 0, 1}; }

   /**
    * @brief Creates a 2D translation transformation matrix
    * @param tx Translation in x direction
    * @param ty Translation in y direction
    * @return Translation transformation matrix
    */
   static Mat3 translation(scalar tx, scalar ty) { return {1, 0, 0, 0, 1, 0, tx, ty, 1}; }

   /**
    * @brief Creates a 2D rotation transformation matrix
    * @param a Rotation angle in radians
    * @return Rotation transformation matrix
    */
   static Mat3 rotation(scalar a) { return {cosf(a), -sinf(a), 0, sinf(a), cosf(a), 0, 0, 0, 1}; }

   /**
    * @brief Creates a 3D rotation matrix around an arbitrary axis
    * @param axis Rotation axis vector (will be normalized)
    * @param angle Rotation angle in radians
    * @return Rotation matrix around the specified axis
    */
   static Mat3 rotationAroundAxis(Vec3 axis, scalar angle) {
      axis = axis.normalized();

      scalar c = cosf(angle);
      scalar s = sinf(angle);
      scalar c_c = 1 - c;
      scalar s_c = 1 - s;

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
   ///@}

   /**
    * @name Accessors
    * Methods for accessing and modifying matrix elements
    */
   ///@{

   /**
    * @brief Gets a column vector from the matrix
    * @param j Column index (0-2)
    * @return Vector representing the specified column
    */
   Vec3 getCol(size_t j) const { return {d[3 * j], d[3 * j + 1], d[3 * j + 2]}; }

   /**
    * @brief Gets a row vector from the matrix
    * @param i Row index (0-2)
    * @return Vector representing the specified row
    */
   Vec3 getRow(size_t i) const { return {d[i], d[i + 3], d[i + 6]}; }

   /**
    * @brief Sets a column of the matrix
    * @param j Column index (0-2)
    * @param v Vector to set as the column
    */
   void setCol(size_t j, const Vec3 &v) {
      d[3 * j] = v[0];
      d[3 * j + 1] = v[1];
      d[3 * j + 2] = v[2];
   }

   /**
    * @brief Sets a row of the matrix
    * @param i Row index (0-2)
    * @param v Vector to set as the row
    */
   void setRow(size_t i, const Vec3 &v) {
      d[i] = v[0];
      d[i + 3] = v[1];
      d[i + 6] = v[2];
   }
   ///@}

   /**
    * @name Matrix Operations
    * Core matrix operations and utilities
    */
   ///@{

   /**
    * @brief Gets the total number of elements in the matrix
    * @return Number of elements (always 9 for a 3x3 matrix)
    */
   static size_t size() { return 9; }

   /**
    * @brief Calculates the determinant of the matrix using Sarrus' rule
    * @return Determinant value
    */
   scalar det() {
      return d[0] * d[4] * d[8] + d[2] * d[3] * d[7] + d[1] * d[5] * d[6] - d[2] * d[4] * d[6] -
             d[1] * d[3] * d[8] - d[0] * d[5] * d[7];
   }

   /**
    * @brief Creates a new matrix that is the transpose of this matrix
    * @return Transposed matrix
    */
   Mat3 transposed() const { return {d[0], d[3], d[6], d[1], d[4], d[7], d[2], d[5], d[8]}; };

   /**
    * @brief Transforms a vector of 3D points by this matrix
    * @param points Vector of 3D points to transform
    * @return Vector of transformed 3D points
    */
   std::vector<Vec3> transformPoints(const std::vector<Vec3> &points) const {
      std::vector<Vec3> transformed_points;
      transformed_points.reserve(points.size());
      for (auto p : points) transformed_points.push_back(*this * p);
      return transformed_points;
   }

   /**
    * @brief Transforms a vector of 2D points by this matrix as a homogeneous transformation
    * @param points Vector of 2D points to transform
    * @return Vector of transformed 2D points
    */
   std::vector<Vec2> transformPoints(const std::vector<Vec2> &points) const {
      std::vector<Vec2> transformed_points;
      transformed_points.reserve(points.size());
      for (auto p : points) transformed_points.push_back(*this * p);
      return transformed_points;
   }

   /**
    * @brief Prints the matrix to stdout in a readable format
    */
   void print() const {
      printf("Mat3:\n");
      for (size_t i = 0; i < 3; i++)
         printf("[%.4f][%.4f][%.4f]\n", d[3 * i], d[3 * i + 1], d[3 * i + 2]);
   }
   ///@}

   /**
    * @name Operators
    * Operator overloads for matrix operations
    */
   ///@{

   /**
    * @brief Array subscript operator (const version)
    * @param i Index of element (0-8)
    * @return Value at the specified index
    */
   scalar operator[](size_t i) const { return d[i]; }

   /**
    * @brief Array subscript operator
    * @param i Index of element (0-8)
    * @return Reference to the element at the specified index
    */
   scalar &operator[](size_t i) { return d[i]; }

   /**
    * @brief Matrix element access operator (const version)
    * @param i Row index (0-2)
    * @param j Column index (0-2)
    * @return Value at the specified position
    */
   scalar operator()(size_t i, size_t j) const { return d[i * 3 + j]; }

   /**
    * @brief Matrix element access operator
    * @param i Row index (0-2)
    * @param j Column index (0-2)
    * @return Reference to the element at the specified position
    */
   scalar &operator()(size_t i, size_t j) { return d[i * 3 + j]; }

   /**
    * @brief Unary minus operator
    * @return Negated matrix
    */
   Mat3 operator-() const {
      return {-d[0], -d[1], -d[2], -d[3], -d[4], -d[5], -d[6], -d[7], -d[8]};
   }

   /**
    * @brief Scalar addition operator
    * @param k Scalar to add to each element
    * @return Result of the addition
    */
   Mat3 operator+(const scalar k) const {
      return {
          d[0] + k, d[1] + k, d[2] + k, d[3] + k, d[4] + k, d[5] + k, d[6] + k, d[7] + k, d[8] + k};
   }

   /**
    * @brief Scalar subtraction operator
    * @param k Scalar to subtract from each element
    * @return Result of the subtraction
    */
   Mat3 operator-(const scalar k) const {
      return {
          d[0] - k, d[1] - k, d[2] - k, d[3] - k, d[4] - k, d[5] - k, d[6] - k, d[7] - k, d[8] - k};
   }

   /**
    * @brief Scalar addition assignment operator
    * @param k Scalar to add to each element
    * @return Reference to this matrix
    */
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

   /**
    * @brief Scalar subtraction assignment operator
    * @param k Scalar to subtract from each element
    * @return Reference to this matrix
    */
   Mat3 &operator-=(const scalar k) {
      d[0] -= k;
      d[1] -= k;
      d[2] -= k;
      d[3] -= k;
      d[4] -= k;
      d[5] -= k;
      d[6] -= k;
      d[7] -= k;
      d[8] -= k;
      return *this;
   }

   /**
    * @brief Scalar multiplication operator
    * @param k Scalar to multiply with each element
    * @return Result of the multiplication
    */
   Mat3 operator*(const scalar k) const {
      return {
          d[0] * k, d[1] * k, d[2] * k, d[3] * k, d[4] * k, d[5] * k, d[6] * k, d[7] * k, d[8] * k};
   }

   /**
    * @brief Scalar division operator
    * @param k Scalar to divide each element by
    * @return Result of the division
    */
   Mat3 operator/(const scalar k) const { return *this * (1 / k); }

   /**
    * @brief Scalar multiplication assignment operator
    * @param k Scalar to multiply with each element
    * @return Reference to this matrix
    */
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

   /**
    * @brief Scalar division assignment operator
    * @param k Scalar to divide each element by
    * @return Reference to this matrix
    */
   Mat3 &operator/=(const scalar k) {
      *this *= (1 / k);
      return *this;
   }

   /**
    * @brief Vector multiplication operator for 3D vectors
    * @param v Vector to multiply with this matrix
    * @return Resulting transformed 3D vector
    */
   Vec3 operator*(const Vec3 &v) const {
      return {d[0] * v.x() + d[3] * v.y() + d[6] * v.z(),
              d[1] * v.x() + d[4] * v.y() + d[7] * v.z(),
              d[2] * v.x() + d[5] * v.y() + d[8] * v.z()};
   }

   /**
    * @brief Vector multiplication operator for 2D vectors (affine transformation)
    * @param v Vector to multiply with this matrix
    * @return Resulting transformed 2D vector
    */
   Vec2 operator*(const Vec2 &v) const {
      scalar w = d[2] * v.x() + d[5] * v.y() + d[8];
      if (w == 0) w = 1;  // Prevent division by 0
      return {(d[0] * v.x() + d[3] * v.y() + d[6]) / w, (d[1] * v.x() + d[4] * v.y() + d[7]) / w};
   }

   /**
    * @brief Matrix addition operator
    * @param m Matrix to add to this matrix
    * @return Result of the addition
    */
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

   /**
    * @brief Matrix subtraction operator
    * @param m Matrix to subtract from this matrix
    * @return Result of the subtraction
    */
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

   /**
    * @brief Matrix addition assignment operator
    * @param m Matrix to add to this matrix
    * @return Reference to this matrix
    */
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

   /**
    * @brief Matrix subtraction assignment operator
    * @param m Matrix to subtract from this matrix
    * @return Reference to this matrix
    */
   Mat3 &operator-=(const Mat3 &m) {
      d[0] -= m.d[0];
      d[1] -= m.d[1];
      d[2] -= m.d[2];
      d[3] -= m.d[3];
      d[4] -= m.d[4];
      d[5] -= m.d[5];
      d[6] -= m.d[6];
      d[7] -= m.d[7];
      d[8] -= m.d[8];
      return *this;
   }

   /**
    * @brief Matrix multiplication operator
    * @param m Matrix to multiply with this matrix
    * @return Result of the multiplication
    */
   Mat3 operator*(const Mat3 &m) const {
      return {d[0] * m.d[0] + d[3] * m.d[1] + d[6] * m.d[2],
              d[1] * m.d[0] + d[4] * m.d[1] + d[7] * m.d[2],
              d[2] * m.d[0] + d[5] * m.d[1] + d[8] * m.d[2],
              d[0] * m.d[3] + d[3] * m.d[4] + d[6] * m.d[5],
              d[1] * m.d[3] + d[4] * m.d[4] + d[7] * m.d[5],
              d[2] * m.d[3] + d[5] * m.d[4] + d[8] * m.d[5],
              d[0] * m.d[6] + d[3] * m.d[7] + d[6] * m.d[8],
              d[1] * m.d[6] + d[4] * m.d[7] + d[7] * m.d[8],
              d[2] * m.d[6] + d[5] * m.d[7] + d[8] * m.d[8]};
   }

   /**
    * @brief Matrix multiplication assignment operator
    * @param m Matrix to multiply with this matrix
    * @return Reference to this matrix
    */
   Mat3 &operator*=(const Mat3 &m) {
      // 1st col
      d[0] = d[0] * m.d[0] + d[3] * m.d[1] + d[6] * m.d[2];
      d[1] = d[1] * m.d[0] + d[4] * m.d[1] + d[7] * m.d[2];
      d[2] = d[2] * m.d[0] + d[5] * m.d[1] + d[8] * m.d[2];
      // 2nd col
      d[3] = d[0] * m.d[3] + d[3] * m.d[4] + d[6] * m.d[5];
      d[4] = d[1] * m.d[3] + d[4] * m.d[4] + d[7] * m.d[5];
      d[5] = d[2] * m.d[3] + d[5] * m.d[4] + d[8] * m.d[5];
      // 3rd col
      d[6] = d[0] * m.d[6] + d[3] * m.d[7] + d[6] * m.d[8];
      d[7] = d[1] * m.d[6] + d[4] * m.d[7] + d[7] * m.d[8];
      d[8] = d[2] * m.d[6] + d[5] * m.d[7] + d[8] * m.d[8];

      return *this;
   }

   /**
    * @brief Equality comparison operator
    * @param m Matrix to compare with
    * @return true if matrices are approximately equal, false otherwise
    */
   bool operator==(const Mat3 &m) const {
      return almostEqual(d[0], m.d[0]) && almostEqual(d[1], m.d[1]) && almostEqual(d[2], m.d[2]) &&
             almostEqual(d[3], m.d[3]) && almostEqual(d[4], m.d[4]) && almostEqual(d[5], m.d[5]) &&
             almostEqual(d[6], m.d[6]) && almostEqual(d[7], m.d[7]) && almostEqual(d[8], m.d[8]);
   }

   /**
    * @brief Inequality comparison operator
    * @param m Matrix to compare with
    * @return true if matrices are not approximately equal, false otherwise
    */
   bool operator!=(const Mat3 &m) const {
      return !almostEqual(d[0], m.d[0]) || !almostEqual(d[1], m.d[1]) ||
             !almostEqual(d[2], m.d[2]) || !almostEqual(d[3], m.d[3]) ||
             !almostEqual(d[4], m.d[4]) || !almostEqual(d[5], m.d[5]) ||
             !almostEqual(d[6], m.d[6]) || !almostEqual(d[7], m.d[7]) || !almostEqual(d[8], m.d[8]);
   }
   ///@}
};

}  // namespace amath

#endif /* __MAT3_H__ */
