#pragma once

#include <amath_core.h>
#include <amath_geometry.h>
#include <amath_eq.h>
#include <vector>

namespace amath {

// struct CollisionShape2D {
//    std::vector<Point2> &points;
//    std::vector<int> indexes;

//    CollisionShape2D(std::vector<Point2> &points, std::vector<int> indexes)
//        : points(points), indexes(indexes) {}

//    bool check_collision(const Vec2 &points) {
//       Point2 &p1 = points[0];
//       const int num_points = indexes.size();
//       for (int i = 1; i < num_points - 1; i++) {
//          Point2 &p2 = points[i + 1];
//          // comprobar colision
//          p1 = p2;
//       }

//    }

bool SegmentIntersection(Segment seg1, Segment seg2) {

   EqSystem2 eqsys(Mat2(seg1.v, -seg2.v), seg2.p - seg1.p);
   EqSol2 sol = eqsys.solve();
   if (sol.type == EqType::kIndependent && sol.values.x() >= 0 && sol.values.x() <= 1 &&
       sol.values.y() >= 0 && sol.values.y() <= 1)
      return true;
   else return false;
}

} // namespace amath