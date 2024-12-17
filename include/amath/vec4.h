#ifndef __AMATH_VEC4_H__
#define __AMATH_VEC4_H__

#include <math.h>
#include <stdio.h>
#include "utils.h"

namespace amath {

/**
 * @brief A 4D vector class for homogeneous coordinates and 4D math
 *
 * This class represents a 4D vector and provides various operations including:
 * - Basic arithmetic operations (+, -, *, /)
 * - Vector operations (dot product, cross product, normalization)
 * - Special handling for homogeneous coordinates (w component)
 * - Length calculations
 * - Common direction vectors (e.g., up)
 *
 * The w component is used for homogeneous coordinates:
 * - w = 0: Vector (direction)
 * - w = 1: Point
 * - w ≠ 0: Homogeneous point (needs normalization)
 */
struct Vec4 {
   scalar d[4];  ///< Array storing vector components (x, y, z, w)

   /**
    * @name Constructors
    * Different ways to construct a Vec4 object
    */
   ///@{

   /**
    * @brief Constructs a vector with specified components
    * @param v0 X component
    * @param v1 Y component
    * @param v2 Z component
    * @param v3 W component
    */
   Vec4(scalar v0, scalar v1, scalar v2, scalar v3) : d{v0, v1, v2, v3} {};

   /**
    * @brief Constructs a vector with all components set to the same value
    * @param v Value for all components
    */
   Vec4(scalar v) : d{v, v, v, v} {};

   /**
    * @brief Constructs a Vec4 from a Vec3, setting w component to 0
    * @param v Value for all components
    */
  //  Vec4(Vec3 v) : d{v.x(), v.y(), v.z(), 0} {};

   /**
    * @brief Default constructor - initializes to zero vector
    */
   Vec4() : Vec4(0, 0, 0, 0) {};

   /**
    * @brief Creates a unit vector pointing up (along Y axis)
    * @return Unit vector (0, 1, 0, 0) in vector form (w=0)
    */
   static Vec4 up() { return {0, 1, 0, 0}; }

   /**
    * @brief Creates a vector with NaN components
    * @return Vector with NaN components
    */
   static Vec4 nan() { return {nanf(""), nanf(""), nanf(""), nanf("")}; }
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
    * @brief Gets the w component (const)
    * @return W component value
    */
   scalar w() const { return d[3]; }

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

   /**
    * @brief Gets reference to w component
    * @return Reference to w component
    */
   scalar &w() { return d[3]; }
   ///@}

   /**
    * @name Vector Operations
    * Core vector operations and utilities
    */
   ///@{

   /**
    * @brief Calculates the squared length of the vector
    * @return Squared length of all components including w
    */
   scalar lengthSquared() const { return x() * x() + y() * y() + z() * z() + w() * w(); }

   /**
    * @brief Calculates the length of the vector
    * @return Length of all components including w
    */
   scalar length() const { return sqrtf(lengthSquared()); }

   /**
    * @brief Returns a normalized version of this vector
    *
    * Normalization behavior depends on w component:
    * - For vectors (w=0): Normalizes to unit length, preserving w=0
    * - For points (w≠0): Normalizes to homogeneous point form with w=1
    *
    * @param tolerance Tolerance for zero length check
    * @return Normalized vector or point
    */
   Vec4 normalized(scalar tolerance = AM_EPSILON) const {
      scalar k = 1 / length();
      // If homogeneous vector (w==0), normalize to length = 1
      if (almostZero(d[3])) return {x() * k, y() * k, z() * k, 0};
      // If homogeneous point (w!=1), normalize to w = 1
      return *this * (1.f / d[3]);
   }

   /**
    * @brief Checks if the vector is normalized (has unit length)
    * @param tolerance Tolerance for unit length check
    * @return true if vector is normalized, false otherwise
    */
   bool isNormalized(const scalar tolerance = AM_EPSILON) {
      return almostZero(lengthSquared() - 1, tolerance);
   }

   /**
    * @brief Calculates dot product with another vector
    * @param v Other vector
    * @return 4D dot product result including w component
    */
   scalar dot(const Vec4 &v) const { return x() * v.x() + y() * v.y() + z() * v.z() + w() * v.w(); }

   /**
    * @brief Calculates cross product with another vector
    * @param v Other vector
    * @return Cross product result (only defined for x,y,z components, w=0)
    * @note This treats the vectors as 3D vectors, ignoring w and setting w=0 in result
    */
   Vec4 cross(const Vec4 &v) const {
      return {y() * v.z() - z() * v.y(), z() * v.x() - x() * v.z(), x() * v.y() - y() * v.x(), 0};
   }

   /**
    * @brief Prints the vector to stdout
    */
   void print() const { printf("Vec4 [%.4f, %.4f, %.4f, %.4f]\n", x(), y(), z(), w()); }
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
   Vec4 operator-() const { return {-d[0], -d[1], -d[2], -d[3]}; }

   /**
    * @brief Scalar addition operator
    * @param k Scalar to add to each component
    * @return Result of the addition
    */
   Vec4 operator+(const scalar k) const { return {x() + k, y() + k, z() + k, w() + k}; }

   /**
    * @brief Scalar subtraction operator
    * @param k Scalar to subtract from each component
    * @return Result of the subtraction
    */
   Vec4 operator-(const scalar k) const { return {x() - k, y() - k, z() - k, w() - k}; }

   /**
    * @brief Scalar multiplication operator
    * @param k Scalar to multiply with each component
    * @return Result of the multiplication
    */
   Vec4 operator*(const scalar k) const { return {x() * k, y() * k, z() * k, w() * k}; }

   /**
    * @brief Scalar multiplication assignment operator
    * @param k Scalar to multiply with each component
    * @return Reference to this vector
    */
   Vec4 &operator*=(const scalar k) {
      d[0] *= k;
      d[1] *= k;
      d[2] *= k;
      d[3] *= k;
      return *this;
   }

   /**
    * @brief Scalar division operator
    * @param k Scalar to divide each component by
    * @return Result of the division
    */
   Vec4 operator/(const scalar k) const { return {x() / k, y() / k, z() / k, w() / k}; }

   /**
    * @brief Scalar division assignment operator
    * @param k Scalar to divide each component by
    * @return Reference to this vector
    */
   Vec4 &operator/=(const scalar k) {
      d[0] /= k;
      d[1] /= k;
      d[2] /= k;
      d[3] /= k;
      return *this;
   }

   /**
    * @brief Vector addition operator
    * @param v Vector to add
    * @return Result of the addition
    */
   Vec4 operator+(const Vec4 &v) const {
      return {x() + v.x(), y() + v.y(), z() + v.z(), w() + v.w()};
   }

   /**
    * @brief Vector subtraction operator
    * @param v Vector to subtract
    * @return Result of the subtraction
    */
   Vec4 operator-(const Vec4 &v) const {
      return {x() - v.x(), y() - v.y(), z() - v.z(), w() - v.w()};
   }

   /**
    * @brief Vector addition assignment operator
    * @param v Vector to add
    * @return Reference to this vector
    */
   Vec4 &operator+=(Vec4 &v) {
      x() += v.x();
      y() += v.y();
      z() += v.z();
      w() += v.w();
      return *this;
   }

   /**
    * @brief Vector subtraction assignment operator
    * @param v Vector to subtract
    * @return Reference to this vector
    */
   Vec4 &operator-=(Vec4 &v) {
      x() -= v.x();
      y() -= v.y();
      z() -= v.z();
      w() -= v.w();
      return *this;
   }

   /**
    * @brief Equality comparison operator
    * @param v Vector to compare with
    * @return true if vectors are approximately equal, false otherwise
    */
   bool operator==(const Vec4 &v) const {
      return fabs(x() - v.x()) < AM_EPSILON && fabs(y() - v.y()) < AM_EPSILON &&
             fabs(z() - v.z()) < AM_EPSILON && fabs(w() - v.w()) < AM_EPSILON;
   }

   /**
    * @brief Inequality comparison operator
    * @param v Vector to compare with
    * @return true if vectors are not approximately equal, false otherwise
    */
   bool operator!=(const Vec4 &v) const {
      return fabs(x() - v.x()) > AM_EPSILON || fabs(y() - v.y()) > AM_EPSILON ||
             fabs(z() - v.z()) > AM_EPSILON || fabs(w() - v.w()) > AM_EPSILON;
   }

   /**
    * @brief Array subscript operator (const version)
    * @param i Index of component (0 for x, 1 for y, 2 for z, 3 for w)
    * @return Value of the specified component
    */
   scalar operator[](const size_t i) const { return d[i]; }

   /**
    * @brief Array subscript operator
    * @param i Index of component (0 for x, 1 for y, 2 for z, 3 for w)
    * @return Reference to the specified component
    */
   scalar &operator[](const size_t i) { return d[i]; }
   ///@}

  //  operator amath::Vec3() const { return amath::Vec3(x(), y(), w()); }
};

}  // namespace amath

#endif /* __AMATH_VEC4_H__ */