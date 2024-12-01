#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <math.h>
#include <stdio.h>
#include <vector>

#include "core.h"
#include "utils.h"
#include "equation.h"

namespace amath {

typedef Vec2 Point2;
typedef Vec3 Point3;

struct Line2 {
   Point2 p;
   Vec2 v;

   Line2(Point2 p, Vec2 v) : p(p), v(v) {}

   scalar distance(const Point2 &point) {
      // TODO
   }

   bool intersect(Line2 other) {
      EqSystem2 eqsys(Mat2(v, -other.v), other.p - p);
      EqSol2 sol = eqsys.solve();
      if (sol.type == EqType::kIndependent && sol.values.x() >= 0 && sol.values.x() <= 1 &&
          sol.values.y() >= 0 && sol.values.y() <= 1)
         return true;
      else return false;
   }
};

struct Line {
   Point3 p;
   Vec3 v;

   Line(Point3 p, Vec3 v) : p(p), v(v) {}

   scalar distance(const Point3 &point) {
      // TODO
   }
};

struct Plane {
   scalar A, B, C, D;  // Ax + By + Cz + D = 0

   Plane(Point3 pt0, Point3 pt1, Point3 pt2, bool ccw_winding = false) {
      Vec3 normal = (pt2 - pt0).cross(pt1 - pt0).normalized();
      normal *= (1 - 2 * ccw_winding);  // Reverse normal if ccw winding

      A = normal.x();
      B = normal.y();
      C = normal.z();
      D = -(A * pt0.x() + B * pt0.y() + C * pt0.z());
   }

   Plane(Vec3 normal, scalar distance) : A(normal.x()), B(normal.y()), C(normal.z()), D(distance) {}

   Vec3 normal() const { return {A, B, C}; }

   scalar distance(Point3 pt) {
      return fabs(A * pt.x() + B * pt.y() + C * pt.z() + D) / sqrtf(A * A + B * B + C * C);
   }

   int side(Point3 pt) { return A * pt.x() + B * pt.y() + C * pt.z() + D > 0 ? 1 : -1; }

   Point3 intersection(Line s) {
      // TODO
   }

   void print() { printf("%+.4fx %+.4fy %+.4fz %+.4f = 0\n", A, B, C, D); }
};

struct RegularPolygon {
   Vec2 center;
   int num_vertex;
   scalar scale;
   scalar rotation;
   bool star;
   scalar star_ratio;

   RegularPolygon()
       : center({0, 0}), num_vertex(0), scale(0), rotation(0), star(0), star_ratio(0) {}

   RegularPolygon(Vec2 center, int num_vertex, scalar scale, scalar rotation = 0, bool star = false,
                  scalar star_ratio = 0.5f)
       : center(center),
         num_vertex(num_vertex),
         scale(scale),
         rotation(rotation),
         star(star),
         star_ratio(star_ratio) {};

   std::vector<scalar> make_path() {
      num_vertex = star ? num_vertex * 2 : num_vertex;

      scalar cos_a = cosf(rotation);
      scalar sin_a = sinf(rotation);
      scalar cos_b = cosf(2 * PI / num_vertex);
      scalar sin_b = sinf(2 * PI / num_vertex);

      std::vector<scalar> points(2 * (num_vertex + 1));

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
      printf(
          "Shape:\n- Center: {%.4f, %.4f}\n"
          "- Scale: %.4f\n"
          "- Rotation: %.4f}\n",
          center.x(),
          center.y(),
          scale,
          rotation);
   }
};

}  // namespace amath

#endif /* __GEOMETRY_H__ */
