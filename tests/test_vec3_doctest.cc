#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "mat3.h"
#include "vec3.h"

TEST_CASE("Matrix determinant computation (Mat3)") {
    amath::Mat3 mat1(1, 2, 3, 4, 5, 6, 7, 8, 9); // Determinant should be 0 (singular matrix)
    CHECK(mat1.det() == doctest::Approx(0));

    amath::Mat3 mat2(6, 1, 1, 4, -2, 5, 2, 8, 7); // Example with determinant -306
    CHECK(mat2.det() == doctest::Approx(-306));
}

TEST_CASE("Matrix transposition (Mat3)") {
    amath::Mat3 mat(1, 2, 3, 4, 5, 6, 7, 8, 9);   // Original matrix
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
    amath::Mat3 mat(1, 0, 0, 0, 1, 0, 0, 0, 1);   // Identity matrix
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

    CHECK(result == doctest::Approx(12)); // 1*4 + 2*(-5) + 3*6 = 12
}

TEST_CASE("Cross product computation (Vec3)") {
    amath::Vec3 vec1(1, 0, 0);
    amath::Vec3 vec2(0, 1, 0);
    amath::Vec3 result = vec1.cross(vec2);

    CHECK(result.x() == 0);
    CHECK(result.y() == 0);
    CHECK(result.z() == 1);
}