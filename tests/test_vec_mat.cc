#include "doctest.h"
#include "amath/linalg.h"

TEST_CASE("Matrix determinant computation") {
   amath::Mat2 mat1(1, 2, 3, 4);  // [1, 3]
                                  // [2, 4]
   CHECK(mat1.det() == doctest::Approx(-2));

   amath::Mat2 mat2(2, 0, 0, 3);  // [2, 0]
                                  // [0, 3]
   CHECK(mat2.det() == doctest::Approx(6));

   amath::Mat2 mat3(0, 0, 0, 0);  // Zero matrix
   CHECK(mat3.det() == doctest::Approx(0));
}

TEST_CASE("Matrix transposition") {
   amath::Mat2 mat(1, 2, 3, 4);  // [1, 3]
                                 // [2, 4]
   amath::Mat2 transposed = mat.transposed();

   CHECK(transposed[0] == 1);
   CHECK(transposed[1] == 3);
   CHECK(transposed[2] == 2);
   CHECK(transposed[3] == 4);
}

TEST_CASE("Matrix-vector multiplication") {
   amath::Mat2 mat(1, 2, 3, 4);  // [1, 3]
                                 // [2, 4]

   amath::Vec2 vec(1, 1);  // [1]
                           // [1]

   amath::Vec2 result = mat * vec;  // [4]
                                    // [6]
   CHECK(result.x() == doctest::Approx(4));
   CHECK(result.y() == doctest::Approx(6));
}

TEST_CASE("Matrix addition") {
   amath::Mat2 mat1(1, 2, 3, 4);  // [1, 3]
                                  // [2, 4]

   amath::Mat2 mat2(4, 3, 2, 1);  // [4, 2]
                                  // [3, 1]

   amath::Mat2 result = mat1 + mat2;  // [5, 5]
                                      // [5, 5]

   CHECK(result[0] == 5);
   CHECK(result[1] == 5);
   CHECK(result[2] == 5);
   CHECK(result[3] == 5);
}

TEST_CASE("Vector-scalar division") {
   amath::Vec2 vec(6, 8);
   amath::Vec2 result = vec / 2;

   CHECK(result.x() == 3);
   CHECK(result.y() == 4);
}

TEST_CASE("Dot product computation") {
   amath::Vec2 vec1(1, 2);
   amath::Vec2 vec2(3, 4);
   double result = vec1.dot(vec2);

   CHECK(result == doctest::Approx(11));
}

TEST_CASE("Matrix determinant computation (Mat3)") {
   amath::Mat3 mat1(1, 2, 3, 4, 5, 6, 7, 8, 9);  // Determinant should be 0 (singular matrix)
   CHECK(mat1.det() == doctest::Approx(0));

   amath::Mat3 mat2(6, 1, 1, 4, -2, 5, 2, 8, 7);  // Example with determinant -306
   CHECK(mat2.det() == doctest::Approx(-306));
}

TEST_CASE("Matrix transposition (Mat3)") {
   amath::Mat3 mat(1, 2, 3, 4, 5, 6, 7, 8, 9);  // Original matrix
   amath::Mat3 transposed = mat.transposed();

   CHECK(transposed[0] == 1);
   CHECK(transposed[1] == 4);
   CHECK(transposed[2] == 7);
   CHECK(transposed[3] == 2);
   CHECK(transposed[4] == 5);
   CHECK(transposed[5] == 8);
   CHECK(transposed[6] == 3);
   CHECK(transposed[7] == 6);
   CHECK(transposed[8] == 9);
}

TEST_CASE("Matrix-vector multiplication (Mat3 and Vec3)") {
   amath::Mat3 mat(1, 0, 0, 0, 1, 0, 0, 0, 1);  // Identity matrix
   amath::Vec3 vec(1, 2, 3);
   amath::Vec3 result = mat * vec;

   CHECK(result.x() == doctest::Approx(1));
   CHECK(result.y() == doctest::Approx(2));
   CHECK(result.z() == doctest::Approx(3));

   amath::Mat3 mat2(2, 0, 0, 0, 3, 0, 0, 0, 4);  // Scaling matrix
   result = mat2 * vec;

   CHECK(result.x() == doctest::Approx(2));
   CHECK(result.y() == doctest::Approx(6));
   CHECK(result.z() == doctest::Approx(12));
}

TEST_CASE("Matrix addition (Mat3)") {
   amath::Mat3 mat1(1, 2, 3, 4, 5, 6, 7, 8, 9);
   amath::Mat3 mat2(9, 8, 7, 6, 5, 4, 3, 2, 1);
   amath::Mat3 result = mat1 + mat2;

   CHECK(result[0] == 10);
   CHECK(result[1] == 10);
   CHECK(result[2] == 10);
   CHECK(result[3] == 10);
   CHECK(result[4] == 10);
   CHECK(result[5] == 10);
   CHECK(result[6] == 10);
   CHECK(result[7] == 10);
   CHECK(result[8] == 10);
}

TEST_CASE("Vector addition (Vec3)") {
   amath::Vec3 vec1(1, 2, 3);
   amath::Vec3 vec2(4, 5, 6);
   amath::Vec3 result = vec1 + vec2;

   CHECK(result.x() == 5);
   CHECK(result.y() == 7);
   CHECK(result.z() == 9);
}

TEST_CASE("Vector-scalar multiplication (Vec3)") {
   amath::Vec3 vec(1, 2, 3);
   amath::Vec3 result = vec * 2;

   CHECK(result.x() == 2);
   CHECK(result.y() == 4);
   CHECK(result.z() == 6);
}

TEST_CASE("Dot product computation (Vec3)") {
   amath::Vec3 vec1(1, 2, 3);
   amath::Vec3 vec2(4, -5, 6);
   double result = vec1.dot(vec2);

   CHECK(result == doctest::Approx(12));  // 1*4 + 2*(-5) + 3*6 = 12
}

TEST_CASE("Cross product computation (Vec3)") {
   amath::Vec3 vec1(1, 0, 0);
   amath::Vec3 vec2(0, 1, 0);
   amath::Vec3 result = vec1.cross(vec2);

   CHECK(result.x() == 0);
   CHECK(result.y() == 0);
   CHECK(result.z() == 1);
}

TEST_CASE("Test Mat4 Constructors and Identity") {
   // Test Mat4 default constructor (should initialize to 0)
   amath::Mat4 mat_default;
   for (int i = 0; i < 16; ++i) { CHECK(mat_default[i] == 0); }

   // Test Mat4 identity constructor
   amath::Mat4 mat_identity = amath::Mat4::identity();
   for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) CHECK(mat_identity(i, j) == (i == j ? 1 : 0));
}

TEST_CASE("Test Mat4 Multiplication") {
   amath::Mat4 mat1 = amath::Mat4::identity();
   amath::Mat4 mat2 = amath::Mat4::scaling(2, 2, 2);
   amath::Mat4 result = mat1 * mat2;

   // The result of multiplying an identity matrix by a scaling matrix should be the scaling matrix
   // itself
   CHECK(result(0, 0) == 2);
   CHECK(result(1, 1) == 2);
   CHECK(result(2, 2) == 2);
   CHECK(result(3, 3) == 1);
}

TEST_CASE("Test Mat4 Rotation") {
   amath::scalar angle = amath::PI / 2;
   amath::Mat4 mat_rotation = amath::Mat4::rotationZ(angle);

   // Check that after rotating 90 degrees, the (1,0) element should be 1 and (0,1) should be -1
   CHECK(mat_rotation(0, 0) == doctest::Approx(0).epsilon(0.01));
   CHECK(mat_rotation(0, 1) == doctest::Approx(1).epsilon(0.01));
   CHECK(mat_rotation(1, 0) == doctest::Approx(-1).epsilon(0.01));
   CHECK(mat_rotation(1, 1) == doctest::Approx(0).epsilon(0.01));
}

TEST_CASE("Test Vec4 and amath::Mat4 Transformation") {
   amath::Vec4 vec(1, 0, 0, 1);  // A point in homogeneous coordinates
   amath::Mat4 mat_translation = amath::Mat4::translation(2, 3, 4);
   amath::Vec4 result = mat_translation * vec;

   // After translating the point (1, 0, 0) by (2, 3, 4), the result should be (3, 3, 4)
   CHECK(result.x() == 3);
   CHECK(result.y() == 3);
   CHECK(result.z() == 4);
   CHECK(result.w() == 1);
}

TEST_CASE("Test amath::Mat4 Determinant") {
   amath::Mat4 mat =
       amath::Mat4::scaling(2, 2, 2);  // A simple scaling matrix with a determinant of 8
   CHECK(mat.det() == 8);
}

TEST_CASE("Test amath::Mat4 Adjugate") {
   amath::Mat4 mat = amath::Mat4::scaling(2, 2, 2);  // Simple scaling matrix
   amath::Mat3 adj = mat.getAdjugate(0, 0);

   // For a scaling matrix, the adjugate of an element is a scalar multiple of a matrix
   CHECK(adj.det() == doctest::Approx(4).epsilon(0.01));
}

TEST_CASE("Test amath::Mat4 Transposition") {
   amath::Mat4 mat = amath::Mat4::scaling(1, 2, 3);
   amath::Mat4 mat_transposed = mat.transposed();

   // Check the transpose of the matrix (swapping rows with columns)
   CHECK(mat_transposed(0, 1) == 0);
   CHECK(mat_transposed(0, 2) == 0);
   CHECK(mat_transposed(1, 0) == 0);
   CHECK(mat_transposed(1, 2) == 0);
   CHECK(mat_transposed(2, 0) == 0);
   CHECK(mat_transposed(2, 1) == 0);
   CHECK(mat_transposed(3, 3) == 1);
}

TEST_CASE("Test amath::Mat4 LookAt") {
   amath::Vec3 camera_pos(0, 0, 5);
   amath::Vec3 camera_dir(0, 0, -1);
   amath::Mat4 mat_view = amath::Mat4::view(camera_pos, camera_dir);

   // Check that the view matrix places the camera at (0, 0, 5) and looks at the origin
   CHECK(mat_view(3, 0) == 0);
   CHECK(mat_view(3, 1) == 0);
   CHECK(mat_view(3, 2) == 5);
}

TEST_CASE("Test Vec4 Cross Product") {
   amath::Vec4 v1(1, 0, 0, 0);
   amath::Vec4 v2(0, 1, 0, 0);
   amath::Vec4 result = v1.cross(v2);

   // Cross product of unit vectors should result in a vector along the third axis
   CHECK(result.x() == 0);
   CHECK(result.y() == 0);
   CHECK(result.z() == 1);
   CHECK(result.w() == 0);
}