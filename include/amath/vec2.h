#ifndef __VEC2_H__
#define __VEC2_H__

#include <math.h>
#include <stdio.h>
#include "utils.h"

namespace amath {

/**
 * @brief A 2D vector class for mathematical operations
 * 
 * This class represents a 2D vector and provides various operations including:
 * - Basic arithmetic operations (+, -, *, /)
 * - Vector operations (dot product, normalization)
 * - Perpendicular vector computation
 * - Length calculations
 */
struct Vec2 {
   scalar d[2]; ///< Array storing vector components (x, y)

   /**
    * @name Constructors
    * Different ways to construct a Vec2 object
    */
   ///@{
   
   /**
    * @brief Constructs a vector with specified components
    * @param v0 X component
    * @param v1 Y component
    */
   Vec2(scalar v0, scalar v1) : d{v0, v1} {};
   
   /**
    * @brief Constructs a vector with both components set to the same value
    * @param v Value for both components
    */
   Vec2(scalar v) : d{v, v} {};
   
   /**
    * @brief Default constructor - initializes to zero vector
    */
   Vec2() : Vec2(0, 0) {};

   /**
    * @brief Creates a vector with NaN components
    * @return Vector with NaN components
    */
   static Vec2 nan() { return {nanf(""), nanf("")}; }
   ///@}

   /**
    * @name Component Access
    * Methods for accessing vector components
    */
   ///@{
   
   /**
    * @brief Gets the x component (const)
    * @return X component value
    */
   scalar x() const { return d[0]; }
   
   /**
    * @brief Gets the y component (const)
    * @return Y component value
    */
   scalar y() const { return d[1]; }
   
   /**
    * @brief Gets reference to x component
    * @return Reference to x component
    */
   scalar &x() { return d[0]; }
   
   /**
    * @brief Gets reference to y component
    * @return Reference to y component
    */
   scalar &y() { return d[1]; }
   ///@}

   /**
    * @name Vector Operations
    * Core vector operations and utilities
    */
   ///@{
   
   /**
    * @brief Calculates the squared length of the vector
    * @return Squared length
    */
   scalar lengthSquared() const { return x() * x() + y() * y(); }
   
   /**
    * @brief Calculates the length of the vector
    * @return Vector length
    */
   scalar length() const { return sqrtf(lengthSquared()); }

   /**
    * @brief Returns a normalized (unit length) version of this vector
    * @param tolerance Tolerance for zero length check
    * @return Normalized vector, or zero vector if length is within tolerance of zero
    */
   Vec2 normalized(scalar tolerance = AM_EPSILON) const {
      scalar l = length();
      if (!almostZero(l, tolerance)) return {x() / l, y() / l};
      else return {0, 0};
   }

   /**
    * @brief Checks if the vector is normalized (has unit length)
    * @param tolerance Tolerance for unit length check
    * @return true if vector is normalized, false otherwise
    */
   bool isNormalized(const scalar tolerance = AM_EPSILON) const {
      return almostZero(lengthSquared() - 1, tolerance);
   }

   /**
    * @brief Returns a vector perpendicular to this one
    * @return Vector rotated 90 degrees counter-clockwise
    */
   Vec2 perpendicular() const { return {-y(), x()}; }

   /**
    * @brief Calculates dot product with another vector
    * @param v Other vector
    * @return Dot product result
    */
   inline scalar dot(const Vec2 &v) const { return x() * v.x() + y() * v.y(); }

   /**
    * @brief Prints the vector to stdout
    */
   void print() const { printf("Vec2 [%.4f, %.4f]\n", x(), y()); }
   ///@}

   /**
    * @name Operators
    * Operator overloads for vector operations
    */
   ///@{
   
   /**
    * @brief Unary minus operator
    * @return Negated vector
    */
   Vec2 operator-() const { return {-d[0], -d[1]}; }

   /**
    * @brief Scalar addition operator
    * @param k Scalar to add to each component
    * @return Result of the addition
    */
   Vec2 operator+(const scalar k) const { return {x() + k, y() + k}; }
   
   /**
    * @brief Scalar subtraction operator
    * @param k Scalar to subtract from each component
    * @return Result of the subtraction
    */
   Vec2 operator-(const scalar k) const { return {x() - k, y() - k}; }
   
   /**
    * @brief Scalar multiplication operator
    * @param k Scalar to multiply with each component
    * @return Result of the multiplication
    */
   Vec2 operator*(const scalar k) const { return {x() * k, y() * k}; }
   
   /**
    * @brief Scalar multiplication assignment operator
    * @param k Scalar to multiply with each component
    * @return Reference to this vector
    */
   Vec2 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      return *this;
   }
   
   /**
    * @brief Scalar division operator
    * @param k Scalar to divide each component by
    * @return Result of the division
    */
   Vec2 operator/(const scalar k) const { return {x() / k, y() / k}; }
   
   /**
    * @brief Scalar division assignment operator
    * @param k Scalar to divide each component by
    * @return Reference to this vector
    */
   Vec2 &operator/=(const scalar k) {
      d[0] /= k;
      d[1] /= k;
      return *this;
   }

   /**
    * @brief Vector addition operator
    * @param v Vector to add
    * @return Result of the addition
    */
   Vec2 operator+(const Vec2 &v) const { return {x() + v.x(), y() + v.y()}; }
   
   /**
    * @brief Vector subtraction operator
    * @param v Vector to subtract
    * @return Result of the subtraction
    */
   Vec2 operator-(const Vec2 &v) const { return {x() - v.x(), y() - v.y()}; }
   
   /**
    * @brief Vector addition assignment operator
    * @param v Vector to add
    * @return Reference to this vector
    */
   Vec2 &operator+=(Vec2 &v) {
      x() += v.x();
      y() += v.y();
      return *this;
   }
   
   /**
    * @brief Vector subtraction assignment operator
    * @param v Vector to subtract
    * @return Reference to this vector
    */
   Vec2 &operator-=(Vec2 &v) {
      x() -= v.x();
      y() -= v.y();
      return *this;
   }

   /**
    * @brief Equality comparison operator
    * @param v Vector to compare with
    * @return true if vectors are approximately equal, false otherwise
    */
   bool operator==(const Vec2 &v) const {
      return fabs(x() - v.x()) < AM_EPSILON && fabs(y() - v.y()) < AM_EPSILON;
   }
   
   /**
    * @brief Inequality comparison operator
    * @param v Vector to compare with
    * @return true if vectors are not approximately equal, false otherwise
    */
   bool operator!=(const Vec2 &v) const {
      return fabs(x() - v.x()) > AM_EPSILON || fabs(y() - v.y()) > AM_EPSILON;
   }

   /**
    * @brief Array subscript operator (const version)
    * @param i Index of component (0 for x, 1 for y)
    * @return Value of the specified component
    */
   scalar operator[](const size_t i) const { return d[i]; }
   
   /**
    * @brief Array subscript operator
    * @param i Index of component (0 for x, 1 for y)
    * @return Reference to the specified component
    */
   scalar &operator[](const size_t i) { return d[i]; }
   ///@}
};

}  // namespace amath

#endif /* __VEC2_H__ */