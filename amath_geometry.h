#pragma once

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <cmath>
#include <cstdio>
#include <vector>

#include <amath_core.h>
#include <amath_utils.h>

namespace amath {

float distance_point_line(const Point2 &point, const Point2 &plane_point,
                          const Vec2 &plane_normal) {
   Vec2 v = plane_point - point;
   return dot_product(v, plane_normal);
}

// float distance_point_plane(const Point3 &point, const Point3 &plane_point,
//                            const Vec3 &plane_normal) {}

Mat3 combine_transforms(std::vector<Mat3> transforms) {
   Mat3 result = Mat3::identity();
   for (auto tr : transforms) result = mat_mul(tr, result);
   return result;
}

struct RegularPolygon {
   Vec2 center;
   int num_vertex;
   float scale;
   float rotation;
   bool star;
   float star_ratio;

   RegularPolygon()
       : center({0, 0}), num_vertex(0), scale(0), rotation(0), star(0), star_ratio(0) {}

   RegularPolygon(Vec2 center, int num_vertex, float scale, float rotation = 0, bool star = false,
                  float star_ratio = 0.5f)
       : center(center), num_vertex(num_vertex), scale(scale), rotation(rotation), star(star),
         star_ratio(star_ratio){};

   std::vector<float> make_path() {
      num_vertex = star ? num_vertex * 2 : num_vertex;

      float cos_a = cosf(rotation);
      float sin_a = sinf(rotation);
      float cos_b = cosf(2 * PI / num_vertex);
      float sin_b = sinf(2 * PI / num_vertex);

      std::vector<float> points(2 * (num_vertex + 1));

      for (int i = 0; i < num_vertex + 1; i++) {
         points[2 * i] = cos_a * cos_b - sin_a * sin_b;
         points[2 * i + 1] = sin_a * cos_b + cos_a * sin_b;

         cos_a = points[2 * i];
         sin_a = points[2 * i + 1];

         points[2 * i] = cos_a * scale + center.x();
         points[2 * i + 1] = sin_a * scale + center.y();
      }

      return points;
   }

   void print() {
      fmt::print("Shape:\n- Center: {{}, {}}\n"
                 "- Scale: {}\n"
                 "- Rotation: {}\n",
                 center.x(),
                 center.y(),
                 scale,
                 rotation);
   }
};

struct Segment {
   Point2 p;
   Vec2 v;

   Segment(Point2 p, Vec2 v) : p(p), v(v) {}
};

struct Plane {
   // Ax + By + Cz + D = 0
   float A, B, C, D;

   Plane(Vec3 pt0, Vec3 pt1, Vec3 pt2) {
      Vec3 normal = cross_product(pt1 - pt0, pt2 - pt0).normalized();
      A = normal.x();
      B = normal.y();
      C = normal.z();
      D = -A * pt0.x() + B * pt0.y() + C * pt0.z();
   }

   Vec3 normal() { return {A, B, C}; }

   void print() { fmt::print("{:+.3f}x {:+.3f}y {:+.3f}z {:+.3f} = 0\n", A, B, C, D); }
};

} // namespace amath
