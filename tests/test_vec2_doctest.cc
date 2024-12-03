#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "mat2.h"
#include "vec2.h"

TEST_CASE("Matrix determinant computation") {
    amath::Mat2 mat1(1, 2, 3, 4); // [1, 3]
                                  // [2, 4]
    CHECK(mat1.det() == doctest::Approx(-2));

    amath::Mat2 mat2(2, 0, 0, 3); // [2, 0]
                                  // [0, 3]
    CHECK(mat2.det() == doctest::Approx(6));

    amath::Mat2 mat3(0, 0, 0, 0); // Zero matrix
    CHECK(mat3.det() == doctest::Approx(0));
}

TEST_CASE("Matrix transposition") {
    amath::Mat2 mat(1, 2, 3, 4);   // [1, 3]
                                   // [2, 4]
    amath::Mat2 transposed = mat.transposed();

    CHECK(transposed[0] == 1);
    CHECK(transposed[1] == 3);
    CHECK(transposed[2] == 2);
    CHECK(transposed[3] == 4);
}

TEST_CASE("Matrix-vector multiplication") {
    amath::Mat2 mat(1, 2, 3, 4);     // [1, 3]
                                     // [2, 4]
                                     
    amath::Vec2 vec(1, 1);           // [1]
                                     // [1]

    amath::Vec2 result = mat * vec;  // [4]
                                     // [6]
    CHECK(result.x() == doctest::Approx(4));
    CHECK(result.y() == doctest::Approx(6));
}

TEST_CASE("Matrix addition") {
    amath::Mat2 mat1(1, 2, 3, 4);    // [1, 3]
                                     // [2, 4]

    amath::Mat2 mat2(4, 3, 2, 1);    // [4, 2]
                                     // [3, 1]

    amath::Mat2 result = mat1 + mat2; // [5, 5]
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
