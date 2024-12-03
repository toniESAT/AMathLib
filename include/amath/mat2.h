#ifndef __MAT2_H__
#define __MAT2_H__

#include <math.h>
#include <stdio.h>
#include <vector>

#include "utils.h"
#include "vec2.h"

/**
 * @brief Mathematical utilities and linear algebra implementations
 */
namespace amath {

/**
 * @brief A 2x2 matrix class for linear algebra operations
 *
 * This class represents a 2x2 matrix and provides various operations including:
 * - Basic arithmetic operations (+, -, *, /)
 * - Matrix-vector multiplication
 * - Matrix-matrix multiplication
 * - Determinant calculation
 * - Transpose operation
 * - Point transformation
 */
struct Mat2 {
   scalar d[4];  ///< Array storing matrix elements in column-major order

   /**
    * @name Constructors
    * Different ways to construct a Mat2 object
    */
   ///@{

   /**
    * @brief Constructs a matrix with specified elements
    * @param m0 Element at position (0,0)
    * @param m1 Element at position (1,0)
    * @param m2 Element at position (0,1)
    * @param m3 Element at position (1,1)
    */
   Mat2(scalar m0, scalar m1, scalar m2, scalar m3) : d{m0, m1, m2, m3} {};

   /**
    * @brief Constructs a matrix with all elements set to the same value
    * @param k Value to set all elements to
    */
   Mat2(scalar k) : Mat2(k, k, k, k) {}

   /**
    * @brief Default constructor - initializes all elements to zero
    */
   Mat2() : d{0} {};

   /**
    * @brief Constructs a matrix from two column vectors
    * @param v0 First column vector
    * @param v1 Second column vector
    */
   Mat2(const Vec2 &v0, const Vec2 &v1) {
      setCol(0, v0);
      setCol(1, v1);
   }
   ///@}

   /**
    * @name Static Factory Methods
    * Methods for creating special matrices
    */
   ///@{

   /**
    * @brief Creates a 2x2 identity matrix
    * @return Identity matrix
    */
   static Mat2 identity() { return {1, 0, 0, 1}; }

   /**
    * @brief Creates a matrix with all elements set to NaN
    * @return Matrix with NaN elements
    */
   static Mat2 nan() { return Mat2(nanf("")); }
   ///@}

   /**
    * @name Accessors
    * Methods for accessing and modifying matrix elements
    */
   ///@{

   /**
    * @brief Gets a column vector from the matrix
    * @param j Column index (0 or 1)
    * @return Vector representing the specified column
    */
   Vec2 getCol(const size_t j) const { return {d[2 * j], d[2 * j + 1]}; }

   /**
    * @brief Gets a row vector from the matrix
    * @param i Row index (0 or 1)
    * @return Vector representing the specified row
    */
   Vec2 getRow(const size_t i) const { return {d[i], d[i + 2]}; }

   /**
    * @brief Sets a column of the matrix
    * @param j Column index (0 or 1)
    * @param v Vector to set as the column
    */
   void setCol(const size_t j, const Vec2 &v) {
      d[2 * j] = v[0];
      d[2 * j + 1] = v[1];
   }

   /**
    * @brief Sets a row of the matrix
    * @param i Row index (0 or 1)
    * @param v Vector to set as the row
    */
   void setRow(const size_t i, const Vec2 &v) {
      d[i] = v[0];
      d[i + 2] = v[1];
   }
   ///@}

   /**
    * @name Matrix Operations
    * Core matrix operations and utilities
    */
   ///@{

   /**
    * @brief Gets the total number of elements in the matrix
    * @return Number of elements (always 4 for a 2x2 matrix)
    */
   static size_t size() { return 4; }

   /**
    * @brief Calculates the determinant of the matrix
    * @return Determinant value
    */
   scalar det() const { return d[0] * d[3] - d[1] * d[2]; }

   /**
    * @brief Creates a new matrix that is the transpose of this matrix
    * @return Transposed matrix
    */
   Mat2 transposed() const { return {d[0], d[2], d[1], d[3]}; };

   /**
    * @brief Transforms a vector of points by this matrix
    * @param points Vector of points to transform
    * @return Vector of transformed points
    */
   std::vector<Vec2> transformPoints(const std::vector<Vec2> &points) const {
      std::vector<Vec2> transformed_points;
      transformed_points.reserve(points.size());
      for (auto p : points) transformed_points.push_back(*this * p);
      return transformed_points;
   }

   /**
    * @brief Prints the matrix to stdout
    */
   void print() const { printf("Mat2:\n[%.4f][%.4f]\n[%.4f][%.4f]", d[0], d[2], d[1], d[3]); }
   ///@}

   /**
    * @name Operators
    * Operator overloads for matrix operations
    */
   ///@{

   /**
    * @brief Array subscript operator (const version)
    * @param i Index of element (0-3)
    * @return Value at the specified index
    */
   scalar operator[](size_t i) const { return d[i]; }

   /**
    * @brief Array subscript operator
    * @param i Index of element (0-3)
    * @return Reference to the element at the specified index
    */
   scalar &operator[](size_t i) { return d[i]; }

   /**
    * @brief Matrix element access operator (const version)
    * @param i Row index (0-1)
    * @param j Column index (0-1)
    * @return Value at the specified position
    */
   scalar operator()(size_t i, size_t j) const {
      if (i >= 0 && i < 2 && j >= 0 && j < 2) return d[i * 2 + j];
   }

   /**
    * @brief Matrix element access operator
    * @param i Row index (0-1)
    * @param j Column index (0-1)
    * @return Reference to the element at the specified position
    */
   scalar &operator()(size_t i, size_t j) {
      if (i >= 0 && i < 2 && j >= 0 && j < 2) return d[i * 2 + j];
   }

   /**
    * @brief Unary minus operator
    * @return Negated matrix
    */
   Mat2 operator-() const { return {-d[0], -d[1], -d[2], -d[3]}; }

   /**
    * @brief Scalar addition operator
    * @param k Scalar to add to each element
    * @return Result of the addition
    */
   Mat2 operator+(const scalar k) const { return {d[0] + k, d[1] + k, d[2] + k, d[3] + k}; }

   /**
    * @brief Scalar subtraction operator
    * @param k Scalar to subtract from each element
    * @return Result of the subtraction
    */
   Mat2 operator-(const scalar k) const { return *this + (-k); }

   /**
    * @brief Scalar addition assignment operator
    * @param k Scalar to add to each element
    * @return Reference to this matrix
    */
   Mat2 &operator+=(const scalar k) {
      d[0] += k;
      d[1] += k;
      d[2] += k;
      d[3] += k;
      return *this;
   }

   /**
    * @brief Scalar subtraction assignment operator
    * @param k Scalar to subtract from each element
    * @return Reference to this matrix
    */
   Mat2 &operator-=(const scalar k) {
      *this += (-k);
      return *this;
   }

   /**
    * @brief Scalar multiplication operator
    * @param k Scalar to multiply with each element
    * @return Result of the multiplication
    */
   Mat2 operator*(const scalar k) const { return {d[0] * k, d[1] * k, d[2] * k, d[3] * k}; }

   /**
    * @brief Scalar division operator
    * @param k Scalar to divide each element by
    * @return Result of the division
    */
   Mat2 operator/(const scalar k) const { return *this * (1 / k); }

   /**
    * @brief Scalar multiplication assignment operator
    * @param k Scalar to multiply with each element
    * @return Reference to this matrix
    */
   Mat2 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      return *this;
   }

   /**
    * @brief Scalar division assignment operator
    * @param k Scalar to divide each element by
    * @return Reference to this matrix
    */
   Mat2 &operator/=(const scalar k) {
      *this *= (1 / k);
      return *this;
   }

   /**
    * @brief Vector multiplication operator
    * @param v Vector to multiply with this matrix
    * @return Resulting transformed vector
    */
   Vec2 operator*(const Vec2 &v) const {
      return {d[0] * v.x() + d[2] * v.y(), d[1] * v.x() + d[3] * v.y()};
   }

   /**
    * @brief Matrix addition operator
    * @param m Matrix to add to this matrix
    * @return Result of the addition
    */
   Mat2 operator+(const Mat2 &m) const {
      return {d[0] + m.d[0], d[1] + m.d[1], d[2] + m.d[2], d[3] + m.d[3]};
   }

   /**
    * @brief Matrix subtraction operator
    * @param m Matrix to subtract from this matrix
    * @return Result of the subtraction
    */
   Mat2 operator-(const Mat2 &m) const {
      return {d[0] - m.d[0], d[1] - m.d[1], d[2] - m.d[2], d[3] - m.d[3]};
   }

   /**
    * @brief Matrix addition assignment operator
    * @param m Matrix to add to this matrix
    * @return Reference to this matrix
    */
   Mat2 &operator+=(const Mat2 &m) {
      d[0] += m.d[0];
      d[1] += m.d[1];
      d[2] += m.d[2];
      d[3] += m.d[3];
      return *this;
   }

   /**
    * @brief Matrix subtraction assignment operator
    * @param m Matrix to subtract from this matrix
    * @return Reference to this matrix
    */
   Mat2 &operator-=(const Mat2 &m) {
      d[0] -= m.d[0];
      d[1] -= m.d[1];
      d[2] -= m.d[2];
      d[3] -= m.d[3];
      return *this;
   }

   /**
    * @brief Matrix multiplication operator
    * @param m Matrix to multiply with this matrix
    * @return Result of the multiplication
    */
   Mat2 operator*(const Mat2 &m) {
      return {d[0] * m.d[0] + d[2] * m.d[1],
              d[1] * m.d[0] + d[3] * m.d[1],
              d[0] * m.d[2] + d[2] * m.d[3],
              d[1] * m.d[2] + d[3] * m.d[3]};
   }

   /**
    * @brief Matrix multiplication assignment operator
    * @param m Matrix to multiply with this matrix
    * @return Reference to this matrix
    */
   Mat2 &operator*=(const Mat2 &m) {
      d[0] = d[0] * m.d[0] + d[2] * m.d[1];
      d[1] = d[1] * m.d[0] + d[3] * m.d[1];
      d[2] = d[0] * m.d[2] + d[2] * m.d[3];
      d[3] = d[1] * m.d[2] + d[3] * m.d[3];

      return *this;
   }

   /**
    * @brief Equality comparison operator
    * @param m Matrix to compare with
    * @return true if matrices are approximately equal, false otherwise
    */
   bool operator==(const Mat2 &m) const {
      return almostEqual(d[0], m.d[0]) && almostEqual(d[1], m.d[1]) && almostEqual(d[2], m.d[2]) &&
             almostEqual(d[3], m.d[3]);
   }

   /**
    * @brief Inequality comparison operator
    * @param m Matrix to compare with
    * @return true if matrices are not approximately equal, false otherwise
    */
   bool operator!=(const Mat2 &m) const {
      return !almostEqual(d[0], m.d[0]) || !almostEqual(d[1], m.d[1]) ||
             !almostEqual(d[2], m.d[2]) || !almostEqual(d[3], m.d[3]);
   }
   ///@}
};

}  // namespace amath

#endif /* __MAT2_H__ */