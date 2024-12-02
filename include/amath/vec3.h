#ifndef __VEC3_H__
#define __VEC3_H__

#include <math.h>
#include <stdio.h>
#include "utils.h"

namespace amath {

/**
 * @brief A 3D vector class for mathematical operations
 * 
 * This class represents a 3D vector and provides various operations including:
 * - Basic arithmetic operations (+, -, *, /)
 * - Vector operations (dot product, cross product, normalization)
 * - Length calculations
 * - Common direction vectors (e.g., up)
 */
struct Vec3 {
   scalar d[3]; ///< Array storing vector components (x, y, z)

   /**
    * @name Constructors
    * Different ways to construct a Vec3 object
    */
   ///@{
   
   /**
    * @brief Constructs a vector with specified components
    * @param v0 X component
    * @param v1 Y component
    * @param v2 Z component
    */
   Vec3(scalar v0, scalar v1, scalar v2) : d{v0, v1, v2} {};
   
   /**
    * @brief Constructs a vector with all components set to the same value
    * @param v Value for all components
    */
   Vec3(scalar v) : d{v, v, v} {};
   
   /**
    * @brief Default constructor - initializes to zero vector
    */
   Vec3() : Vec3(0, 0, 0) {};

   /**
    * @brief Creates a unit vector pointing up (along Y axis)
    * @return Unit vector (0, 1, 0)
    */
   static Vec3 up() { return {0, 1, 0}; }
   
   /**
    * @brief Creates a vector with NaN components
    * @return Vector with NaN components
    */
   static Vec3 nan() { return {nanf(""), nanf(""), nanf("")}; }
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
    * @brief Gets the z component (const)
    * @return Z component value
    */
   scalar z() const { return d[2]; }
   
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
   
   /**
    * @brief Gets reference to z component
    * @return Reference to z component
    */
   scalar &z() { return d[2]; }
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
   scalar lengthSquared() const { return x() * x() + y() * y() + z() * z(); }
   
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
   Vec3 normalized(scalar tolerance = AM_EPSILON) const {
      scalar l = length();
      if (!almostZero(l, tolerance)) return {x() / l, y() / l, z() / l};
      else return {0, 0, 0};
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
    * @brief Calculates dot product with another vector
    * @param v Other vector
    * @return Dot product result
    */
   scalar dot(const Vec3 &v) const { return x() * v.x() + y() * v.y() + z() * v.z(); }

   /**
    * @brief Calculates cross product with another vector
    * @param v Other vector
    * @return Cross product result
    */
   Vec3 cross(const Vec3 &v) const {
      return {y() * v.z() - z() * v.y(), z() * v.x() - x() * v.z(), x() * v.y() - y() * v.x()};
   }

   /**
    * @brief Prints the vector to stdout
    */
   void print() const { printf("Vec3 [%.4f,%.4f, %.4f]\n", x(), y(), z()); }
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
   Vec3 operator-() const { return {-d[0], -d[1], -d[2]}; }

   /**
    * @brief Scalar addition operator
    * @param k Scalar to add to each component
    * @return Result of the addition
    */
   Vec3 operator+(const scalar k) const { return {x() + k, y() + k, z() + k}; }
   
   /**
    * @brief Scalar subtraction operator
    * @param k Scalar to subtract from each component
    * @return Result of the subtraction
    */
   Vec3 operator-(const scalar k) const { return {x() - k, y() - k, z() - k}; }
   
   /**
    * @brief Scalar multiplication operator
    * @param k Scalar to multiply with each component
    * @return Result of the multiplication
    */
   Vec3 operator*(const scalar k) const { return {x() * k, y() * k, z() * k}; }
   
   /**
    * @brief Scalar multiplication assignment operator
    * @param k Scalar to multiply with each component
    * @return Reference to this vector
    */
   Vec3 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      return *this;
   }
   
   /**
    * @brief Scalar division operator
    * @param k Scalar to divide each component by
    * @return Result of the division
    */
   Vec3 operator/(const scalar k) const { return {x() / k, y() / k, z() / k}; }
   
   /**
    * @brief Scalar division assignment operator
    * @param k Scalar to divide each component by
    * @return Reference to this vector
    */
   Vec3 &operator/=(const scalar k) {
      d[0] /= k;
      d[1] /= k;
      d[2] /= k;
      return *this;
   }

   /**
    * @brief Vector addition operator
    * @param v Vector to add
    * @return Result of the addition
    */
   Vec3 operator+(const Vec3 &v) const { return {x() + v.x(), y() + v.y(), z() + v.z()}; }
   
   /**
    * @brief Vector subtraction operator
    * @param v Vector to subtract
    * @return Result of the subtraction
    */
   Vec3 operator-(const Vec3 &v) const { return {x() - v.x(), y() - v.y(), z() - v.z()}; }
   
   /**
    * @brief Vector addition assignment operator
    * @param v Vector to add
    * @return Reference to this vector
    */
   Vec3 &operator+=(Vec3 &v) {
      x() += v.x();
      y() += v.y();
      z() += v.z();
      return *this;
   }
   
   /**
    * @brief Vector subtraction assignment operator
    * @param v Vector to subtract
    * @return Reference to this vector
    */
   Vec3 &operator-=(Vec3 &v) {
      x() -= v.x();
      y() -= v.y();
      z() -= v.z(); 
      return *this;
   }

   /**
    * @brief Equality comparison operator
    * @param v Vector to compare with
    * @return true if vectors are approximately equal, false otherwise
    */
   bool operator==(const Vec3 &v) const {
      return fabs(x() - v.x()) < AM_EPSILON && fabs(y() - v.y()) < AM_EPSILON &&
             fabs(z() - v.z()) < AM_EPSILON;
   }
   
   /**
    * @brief Inequality comparison operator
    * @param v Vector to compare with
    * @return true if vectors are not approximately equal, false otherwise
    */
   bool operator!=(const Vec3 &v) const {
      return fabs(x() - v.x()) > AM_EPSILON || fabs(y() - v.y()) > AM_EPSILON ||
             fabs(z() - v.z()) > AM_EPSILON;
   }

   /**
    * @brief Array subscript operator (const version)
    * @param i Index of component (0 for x, 1 for y, 2 for z)
    * @return Value of the specified component
    */
   scalar operator[](const size_t i) const { return d[i]; }
   
   /**
    * @brief Array subscript operator
    * @param i Index of component (0 for x, 1 for y, 2 for z)
    * @return Reference to the specified component
    */
   scalar &operator[](const size_t i) { return d[i]; }
   ///@}
};

}  // namespace amath

#endif /* __VEC3_H__ */