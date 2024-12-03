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

/**
 * @brief A 4x4 matrix class for 3D transformations and projections
 *
 * This class represents a 4x4 matrix and provides various operations including:
 * - Basic arithmetic operations (+, -, *, /)
 * - Matrix-vector multiplication
 * - Matrix-matrix multiplication
 * - Determinant calculation
 * - Transpose operation
 * - Point transformation
 * - Common 3D transformation matrices (rotation, translation, scaling)
 * - View and projection matrices for 3D graphics
 */
struct Mat4 {
   scalar d[16];  ///< Array storing matrix elements in column-major order

   /**
    * @name Constructors
    * Different ways to construct a Mat4 object
    */
   ///@{

   /**
    * @brief Constructs a matrix with specified elements
    * @param m0-m15 Elements in column-major order
    */
   Mat4(scalar m0, scalar m1, scalar m2, scalar m3, scalar m4, scalar m5, scalar m6, scalar m7,
        scalar m8, scalar m9, scalar m10, scalar m11, scalar m12, scalar m13, scalar m14,
        scalar m15)
       : d{m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15} {};

   /**
    * @brief Constructs a matrix with all elements set to the same value
    * @param k Value to set all elements to
    */
   Mat4(scalar k) : Mat4(k, k, k, k, k, k, k, k, k, k, k, k, k, k, k, k) {}

   /**
    * @brief Default constructor - initializes all elements to zero
    */
   Mat4() : d{0} {};

   /**
    * @brief Constructs a matrix from four column vectors
    * @param v0 First column vector
    * @param v1 Second column vector
    * @param v2 Third column vector
    * @param v3 Fourth column vector
    */
   Mat4(const Vec4 &v0, const Vec4 &v1, const Vec4 &v2, const Vec4 &v3) {
      setCol(0, v0);
      setCol(1, v1);
      setCol(2, v2);
      setCol(3, v3);
   }
   ///@}

   /**
    * @name Static Factory Methods - Basic Matrices
    * Methods for creating special matrices
    */
   ///@{

   /**
    * @brief Creates a 4x4 identity matrix
    * @return Identity matrix
    */
   static Mat4 identity() { return {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}; }

   /**
    * @brief Creates a matrix with all elements set to NaN
    * @return Matrix with NaN elements
    */
   static Mat4 nan() { return Mat4(nanf("")); };

   /**
    * @brief Creates a 3D scaling transformation matrix
    * @param sx Scale factor in x direction
    * @param sy Scale factor in y direction
    * @param sz Scale factor in z direction
    * @return Scaling transformation matrix
    */
   static Mat4 scaling(scalar sx, scalar sy, scalar sz) {
      return {sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0, 0, 0, 0, 1};
   }

   /**
    * @brief Creates a 3D translation transformation matrix
    * @param tx Translation in x direction
    * @param ty Translation in y direction
    * @param tz Translation in z direction
    * @return Translation transformation matrix
    */
   static Mat4 translation(scalar tx, scalar ty, scalar tz) {
      return {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, tx, ty, tz, 1};
   }
   ///@}

   /**
    * @name Static Factory Methods - Rotation Matrices
    * Methods for creating various types of rotation matrices
    */
   ///@{

   /**
    * @brief Creates a rotation matrix around the X axis
    * @param a Rotation angle in radians
    * @return Rotation matrix
    */
   static Mat4 rotationX(scalar a) {
      return {1, 0, 0, 0, 0, cosf(a), -sinf(a), 0, 0, sinf(a), cosf(a), 0, 0, 0, 0, 1};
   }

   /**
    * @brief Creates a rotation matrix around the Y axis
    * @param a Rotation angle in radians
    * @return Rotation matrix
    */
   static Mat4 rotationY(scalar a) {
      return {cosf(a), 0, sinf(a), 0, 0, 1, 0, 0, -sinf(a), 0, cosf(a), 0, 0, 0, 0, 1};
   }

   /**
    * @brief Creates a rotation matrix around the Z axis
    * @param a Rotation angle in radians
    * @return Rotation matrix
    */
   static Mat4 rotationZ(scalar a) {
      return {cosf(a), sinf(a), 0, 0, -sinf(a), cosf(a), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
   }

   /**
    * @brief Creates a rotation matrix from Euler angles (XYZ order)
    * @param x Rotation around X axis in radians
    * @param y Rotation around Y axis in radians
    * @param z Rotation around Z axis in radians
    * @return Combined rotation matrix
    */
   static Mat4 rotation(scalar x, scalar y, scalar z) {
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

   /**
    * @brief Creates a rotation matrix around an arbitrary axis
    * @param axis Rotation axis vector (will be normalized)
    * @param angle Rotation angle in radians
    * @return Rotation matrix around the specified axis
    */
   static Mat4 rotationAroundAxis(Vec4 axis, scalar angle) {
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
   ///@}
   /**
    * @name Static Factory Methods - Transform and Camera Matrices
    * Methods for creating transformation and view matrices
    */
   ///@{

   /**
    * @brief Creates a combined transformation matrix (translation * rotation * scale)
    * @param translate Translation vector (default: no translation)
    * @param scale Scale vector (default: no scaling)
    * @param rot Rotation vector in radians (default: no rotation)
    * @return Combined transformation matrix
    */
   static Mat4 transform(Vec3 translate = {0, 0, 0}, Vec3 scale = {1, 1, 1}, Vec3 rot = {0, 0, 0}) {
      return Mat4::translation(translate.x(), translate.y(), translate.z()) *
             (Mat4::rotation(rot.x(), rot.y(), rot.z()) *
              Mat4::scaling(scale.x(), scale.y(), scale.z()));
   }

   /**
    * @brief Creates a perspective projection matrix
    * @param fov Field of view angle in radians (default: 90 degrees)
    * @param aspect Aspect ratio (width/height) (default: 1)
    * @param zNear Near clipping plane distance (default: 1)
    * @param zFar Far clipping plane distance (default: 100)
    * @return Perspective projection matrix
    */
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
              -zNear * zFar / (zFar - zNear),
              0};
   }

   /**
    * @brief Creates a view matrix from camera position and direction
    * @param camera_pos Camera position in world space
    * @param camera_dir Camera direction vector
    * @return View transformation matrix
    */
   static Mat4 view(Vec3 camera_pos, Vec3 camera_dir) {
      Vec3 right = Vec3(0, 1, 0).cross(camera_dir);
      Vec3 up = camera_dir.cross(right);

      return {right.x(),
              up.x(),
              camera_dir.x(),
              0,
              right.y(),
              up.y(),
              camera_dir.y(),
              0,
              right.z(),
              up.z(),
              camera_dir.z(),
              0,
              -right.dot(camera_pos),
              -up.dot(camera_pos),
              -camera_dir.dot(camera_pos),
              1};
   }

   /**
    * @brief Creates a look-at view matrix
    * @param position Camera position in world space
    * @param target Point the camera is looking at
    * @param world_up World's up vector (default: +Z axis)
    * @return Look-at view matrix
    */
   static Mat4 lookAt(Vec3 position, Vec3 target, Vec3 world_up = {0, 0, 1}) {
      Vec3 direction = Vec3(0, 1, 0).cross(target - position);
      Vec3 right = world_up.cross(direction);
      Vec3 up = direction.cross(right);

      return {right.x(),
              up.x(),
              direction.x(),
              0,
              right.y(),
              up.y(),
              direction.y(),
              0,
              right.z(),
              up.z(),
              direction.z(),
              0,
              -right.dot(position),
              -up.dot(position),
              -direction.dot(position),
              1};
   }

   /**
    * @brief Creates an orthographic projection matrix
    * @param left Left clipping plane x-coordinate
    * @param right Right clipping plane x-coordinate
    * @param bottom Bottom clipping plane y-coordinate
    * @param top Top clipping plane y-coordinate
    * @param zNear Near clipping plane z-coordinate
    * @param zFar Far clipping plane z-coordinate
    * @return Orthographic projection matrix
    */
   static Mat4 ortho(scalar left, scalar right, scalar bottom, scalar top, scalar zNear,
                     scalar zFar) {
      return {2 / (right - left),
              0,
              0,
              0,
              0,
              2 / (top - bottom),
              0,
              0,
              0,
              0,
              -2 / (zFar - zNear),
              0,
              -(right + left) / (right - left),
              -(top + bottom) / (top - bottom),
              (zFar + zNear) / (zFar - zNear),
              1};
   }
   ///@}
   /**
    * @name Accessors
    * Methods for accessing and modifying matrix elements
    */
   ///@{

   /**
    * @brief Gets a column vector from the matrix
    * @param j Column index (0-3)
    * @return Vector representing the specified column
    */
   Vec4 getCol(size_t j) const { return {d[4 * j], d[4 * j + 1], d[4 * j + 2], d[4 * j + 3]}; }

   /**
    * @brief Gets a row vector from the matrix
    * @param i Row index (0-3)
    * @return Vector representing the specified row
    */
   Vec4 getRow(size_t i) const { return {d[i], d[i + 4], d[i + 8], d[i + 12]}; }

   /**
    * @brief Sets a column of the matrix
    * @param j Column index (0-3)
    * @param v Vector to set as the column
    */
   void setCol(size_t j, const Vec4 &v) {
      d[4 * j] = v[0];
      d[4 * j + 1] = v[1];
      d[4 * j + 2] = v[2];
      d[4 * j + 3] = v[3];
   }

   /**
    * @brief Sets a row of the matrix
    * @param i Row index (0-3)
    * @param v Vector to set as the row
    */
   void setRow(size_t i, const Vec4 &v) {
      d[i] = v[0];
      d[i + 4] = v[1];
      d[i + 8] = v[2];
      d[i + 12] = v[3];
   }
   ///@}

   /**
    * @name Matrix Operations
    * Core matrix operations and utilities
    */
   ///@{

   /**
    * @brief Gets the total number of elements in the matrix
    * @return Number of elements (always 16 for a 4x4 matrix)
    */
   static size_t size() { return 16; }

   /**
    * @brief Prints the matrix to stdout in a readable format
    */
   void print() const {
      printf("Mat4\n");
      for (size_t i = 0; i < 4; i++)
         printf(
             "[%.4f][%.4f][%.4f][%.4f]\n", d[4 * 0 + i], d[4 * 1 + i], d[4 * 2 + i], d[4 * 3 + i]);
      printf("\n");
   }

   /**
    * @brief Gets the adjugate (adjoint) matrix for a specific element
    * @param row Row index of the element
    * @param col Column index of the element
    * @return 3x3 adjugate matrix
    */
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

   /**
    * @brief Calculates the determinant of the matrix using cofactor expansion
    * @return Determinant value
    */
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

   /**
    * @brief Creates a new matrix that is the transpose of this matrix
    * @return Transposed matrix
    */
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
   }

   /**
    * @brief Transforms a vector of 4D points by this matrix
    * @param points Vector of 4D points to transform
    * @return Vector of transformed 4D points
    */
   std::vector<Vec4> transformPoints(const std::vector<Vec4> &points) const {
      std::vector<Vec4> transformed_points;
      transformed_points.reserve(points.size());
      for (auto p : points) transformed_points.push_back(*this * p);
      return transformed_points;
   }

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
   ///@}
   /**
    * @name Operators
    * Operator overloads for matrix operations
    */
   ///@{

   /**
    * @brief Array subscript operator (const version)
    * @param i Index of element (0-15)
    * @return Value at the specified index
    */
   scalar operator[](size_t i) const { return d[i]; }

   /**
    * @brief Array subscript operator
    * @param i Index of element (0-15)
    * @return Reference to the element at the specified index
    */
   scalar &operator[](size_t i) { return d[i]; }

   /**
    * @brief Matrix element access operator (const version)
    * @param i Row index (0-3)
    * @param j Column index (0-3)
    * @return Value at the specified position
    */
   scalar operator()(size_t i, size_t j) const {
      if (i >= 0 && i < 4 && j >= 0 && j < 4) return d[i * 4 + j];
   }

   /**
    * @brief Matrix element access operator
    * @param i Row index (0-3)
    * @param j Column index (0-3)
    * @return Reference to the element at the specified position
    */
   scalar &operator()(size_t i, size_t j) {
      if (i >= 0 && i < 4 && j >= 0 && j < 4) return d[i * 4 + j];
   }

   /**
    * @brief Unary minus operator
    * @return Negated matrix
    */
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

   /**
    * @brief Scalar addition operator
    * @param k Scalar to add to each element
    * @return Result of the addition
    */
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

   /**
    * @brief Scalar subtraction operator
    * @param k Scalar to subtract from each element
    * @return Result of the subtraction
    */
   Mat4 operator-(scalar k) const { return *this + (-k); }

   /**
    * @brief Scalar addition assignment operator
    * @param k Scalar to add to each element
    * @return Reference to this matrix
    */
   Mat4 &operator+=(scalar k) {
      for (int i = 0; i < 16; i++) d[i] += k;
      return *this;
   }

   /**
    * @brief Scalar subtraction assignment operator
    * @param k Scalar to subtract from each element
    * @return Reference to this matrix
    */
   Mat4 &operator-=(scalar k) {
      *this += (-k);
      return *this;
   }

   /**
    * @brief Scalar multiplication operator
    * @param k Scalar to multiply with each element
    * @return Result of the multiplication
    */
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

   /**
    * @brief Scalar division operator
    * @param k Scalar to divide each element by
    * @return Result of the division
    */
   Mat4 operator/(scalar k) const { return *this * (1 / k); }

   /**
    * @brief Scalar multiplication assignment operator
    * @param k Scalar to multiply with each element
    * @return Reference to this matrix
    */
   Mat4 operator*=(scalar k) {
      for (int i = 0; i < 16; i++) d[i] *= k;
      return *this;
   }

   /**
    * @brief Scalar division assignment operator
    * @param k Scalar to divide each element by
    * @return Reference to this matrix
    */
   Mat4 operator/=(scalar k) {
      *this *= (1 / k);
      return *this;
   }

   /**
    * @brief Vector multiplication operator for 4D vectors
    * @param v Vector to multiply with this matrix
    * @return Resulting transformed vector
    */
   Vec4 operator*(const Vec4 &v) const {
      return {d[0] * v.x() + d[4] * v.y() + d[8] * v.z() + d[12] * v.w(),
              d[1] * v.x() + d[5] * v.y() + d[9] * v.z() + d[13] * v.w(),
              d[2] * v.x() + d[6] * v.y() + d[10] * v.z() + d[14] * v.w(),
              d[3] * v.x() + d[7] * v.y() + d[11] * v.z() + d[15] * v.w()};
   }

   /**
    * @brief Vector multiplication operator for 3D vectors (affine transformation)
    * @param v Vector to multiply with this matrix
    * @return Resulting transformed 3D vector
    */
   Vec3 operator*(const Vec3 &v) const {
      scalar w = d[3] * v.x() + d[7] * v.y() + d[11] * v.z() + d[15];
      if (w == 0) w = 1;  // Prevent division by 0
      return {(d[0] * v.x() + d[4] * v.y() + d[8] * v.z() + d[12]) / w,
              (d[1] * v.x() + d[5] * v.y() + d[9] * v.z() + d[13]) / w,
              (d[2] * v.x() + d[6] * v.y() + d[10] * v.z() + d[14]) / w};
   }

   /**
    * @brief Matrix addition operator
    * @param m Matrix to add to this matrix
    * @return Result of the addition
    */
   Mat4 operator+(const Mat4 &m) const {
      Mat4 result;
      for (int i = 0; i < 16; i++) result.d[i] = d[i] + m.d[i];
      return result;
   }

   /**
    * @brief Matrix subtraction operator
    * @param m Matrix to subtract from this matrix
    * @return Result of the subtraction
    */
   Mat4 operator-(const Mat4 &m) const {
      Mat4 result;
      for (int i = 0; i < 16; i++) result.d[i] = d[i] - m.d[i];
      return result;
   }

   /**
    * @brief Matrix addition assignment operator
    * @param m Matrix to add to this matrix
    * @return Reference to this matrix
    */
   Mat4 &operator+=(const Mat4 &m) {
      for (int i = 0; i < 16; i++) d[i] += m.d[i];
      return *this;
   }

   /**
    * @brief Matrix subtraction assignment operator
    * @param m Matrix to subtract from this matrix
    * @return Reference to this matrix
    */
   Mat4 &operator-=(const Mat4 &m) {
      for (int i = 0; i < 16; i++) d[i] -= m.d[i];
      return *this;
   }

   /**
    * @brief Matrix multiplication operator
    * @param m Matrix to multiply with this matrix
    * @return Result of the multiplication
    */
   Mat4 operator*(const Mat4 &m) {
      return {d[0] * m.d[0] + d[4] * m.d[1] + d[8] * m.d[2] + d[12] * m.d[3],
              d[1] * m.d[0] + d[5] * m.d[1] + d[9] * m.d[2] + d[13] * m.d[3],
              d[2] * m.d[0] + d[6] * m.d[1] + d[10] * m.d[2] + d[14] * m.d[3],
              d[3] * m.d[0] + d[7] * m.d[1] + d[11] * m.d[2] + d[15] * m.d[3],
              d[0] * m.d[4] + d[4] * m.d[5] + d[8] * m.d[6] + d[12] * m.d[7],
              d[1] * m.d[4] + d[5] * m.d[5] + d[9] * m.d[6] + d[13] * m.d[7],
              d[2] * m.d[4] + d[6] * m.d[5] + d[10] * m.d[6] + d[14] * m.d[7],
              d[3] * m.d[4] + d[7] * m.d[5] + d[11] * m.d[6] + d[15] * m.d[7],
              d[0] * m.d[8] + d[4] * m.d[9] + d[8] * m.d[10] + d[12] * m.d[11],
              d[1] * m.d[8] + d[5] * m.d[9] + d[9] * m.d[10] + d[13] * m.d[11],
              d[2] * m.d[8] + d[6] * m.d[9] + d[10] * m.d[10] + d[14] * m.d[11],
              d[3] * m.d[8] + d[7] * m.d[9] + d[11] * m.d[10] + d[15] * m.d[11],
              d[0] * m.d[12] + d[4] * m.d[13] + d[8] * m.d[14] + d[12] * m.d[15],
              d[1] * m.d[12] + d[5] * m.d[13] + d[9] * m.d[14] + d[13] * m.d[15],
              d[2] * m.d[12] + d[6] * m.d[13] + d[10] * m.d[14] + d[14] * m.d[15],
              d[3] * m.d[12] + d[7] * m.d[13] + d[11] * m.d[14] + d[15] * m.d[15]};
   }

   /**
    * @brief Matrix multiplication assignment operator
    * @param m Matrix to multiply with this matrix
    * @return Reference to this matrix
    */
   Mat4 &operator*=(const Mat4 &m) {
      d[0] = d[0] * m.d[0] + d[4] * m.d[1] + d[8] * m.d[2] + d[12] * m.d[3];
      d[1] = d[1] * m.d[0] + d[5] * m.d[1] + d[9] * m.d[2] + d[13] * m.d[3];
      d[2] = d[2] * m.d[0] + d[6] * m.d[1] + d[10] * m.d[2] + d[14] * m.d[3];
      d[3] = d[3] * m.d[0] + d[7] * m.d[1] + d[11] * m.d[2] + d[15] * m.d[3];
      d[4] = d[0] * m.d[4] + d[4] * m.d[5] + d[8] * m.d[6] + d[12] * m.d[7];
      d[5] = d[1] * m.d[4] + d[5] * m.d[5] + d[9] * m.d[6] + d[13] * m.d[7];
      d[6] = d[2] * m.d[4] + d[6] * m.d[5] + d[10] * m.d[6] + d[14] * m.d[7];
      d[7] = d[3] * m.d[4] + d[7] * m.d[5] + d[11] * m.d[6] + d[15] * m.d[7];
      d[8] = d[0] * m.d[8] + d[4] * m.d[9] + d[8] * m.d[10] + d[12] * m.d[11];
      d[9] = d[1] * m.d[8] + d[5] * m.d[9] + d[9] * m.d[10] + d[13] * m.d[11];
      d[10] = d[2] * m.d[8] + d[6] * m.d[9] + d[10] * m.d[10] + d[14] * m.d[11];
      d[11] = d[3] * m.d[8] + d[7] * m.d[9] + d[11] * m.d[10] + d[15] * m.d[11];
      d[12] = d[0] * m.d[12] + d[4] * m.d[13] + d[8] * m.d[14] + d[12] * m.d[15];
      d[13] = d[1] * m.d[12] + d[5] * m.d[13] + d[9] * m.d[14] + d[13] * m.d[15];
      d[14] = d[2] * m.d[12] + d[6] * m.d[13] + d[10] * m.d[14] + d[14] * m.d[15];
      d[15] = d[3] * m.d[12] + d[7] * m.d[13] + d[11] * m.d[14] + d[15] * m.d[15];

      return *this;
   }

   /**
    * @brief Equality comparison operator
    * @param m Matrix to compare with
    * @return true if matrices are approximately equal, false otherwise
    * @todo Could be optimized with SIMD operations
    */
   bool operator==(const Mat3 &m) const {
      return almostEqual(d[0], m.d[0]) && almostEqual(d[1], m.d[1]) && almostEqual(d[2], m.d[2]) &&
             almostEqual(d[3], m.d[3]) && almostEqual(d[4], m.d[4]) && almostEqual(d[5], m.d[5]) &&
             almostEqual(d[6], m.d[6]) && almostEqual(d[7], m.d[7]) && almostEqual(d[8], m.d[8]) &&
             almostEqual(d[9], m.d[9]) && almostEqual(d[10], m.d[10]) &&
             almostEqual(d[11], m.d[11]) && almostEqual(d[12], m.d[12]) &&
             almostEqual(d[13], m.d[13]) && almostEqual(d[14], m.d[14]) &&
             almostEqual(d[15], m.d[15]);
   }

   /**
    * @brief Inequality comparison operator
    * @param m Matrix to compare with
    * @return true if matrices are not approximately equal, false otherwise
    * @todo Could be optimized with SIMD operations
    */
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
   ///@}
};

}  // namespace amath

#endif /* __MAT4_H__ */
