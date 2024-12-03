#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "Mat4.h"
#include "vec4.h"


TEST_CASE("Test Mat4 Constructors and Identity") {
    // Test Mat4 default constructor (should initialize to 0)
    amath::Mat4 mat_default;
    for (int i = 0; i < 16; ++i) {
        CHECK(mat_default[i] == 0);
    }

    // Test Mat4 identity constructor
    amath::Mat4 mat_identity = amath::Mat4::identity();
    CHECK(mat_identity(0, 0) == 1);
    CHECK(mat_identity(1, 1) == 1);
    CHECK(mat_identity(2, 2) == 1);
    CHECK(mat_identity(3, 3) == 1);
    for (int i = 0; i < 3; ++i) {
        CHECK(mat_identity(0, i) == 0);
        CHECK(mat_identity(1, i) == 0);
        CHECK(mat_identity(2, i) == 0);
        CHECK(mat_identity(3, i) == 0);
    }
}

TEST_CASE("Test Mat4 Multiplication") {
    amath::Mat4 mat1 = amath::Mat4::identity();
    amath::Mat4 mat2 = amath::Mat4::scaling(2, 2, 2);
    amath::Mat4 result = mat1 * mat2;

    // The result of multiplying an identity matrix by a scaling matrix should be the scaling matrix itself
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
    amath::Mat4 mat = amath::Mat4::scaling(2, 2, 2);  // A simple scaling matrix with a determinant of 8
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