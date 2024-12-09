#ifndef __AMATH_GEOMETRY_H__
#define __AMATH_GEOMETRY_H__

#include <math.h>
#include <stdio.h>
#include <vector>

#include "linalg.h"
#include "utils.h"
#include "equation.h"

namespace amath {

/** @brief Type alias for 2D point using Vec2 */
typedef Vec2 Point2;
/** @brief Type alias for 3D point using Vec3 */
typedef Vec3 Point3;

/**
 * @brief Represents a 2D line segment or infinite line
 *
 * The line is defined by a point and a direction vector.
 * Can be used for both line segments and infinite lines depending on context.
 */
struct Line2 {
   Point2 p;  ///< Point on the line
   Vec2 v;    ///< Direction vector of the line

   /**
    * @brief Constructs a line from a point and direction
    * @param p Point on the line
    * @param v Direction vector
    */
   Line2(Point2 p, Vec2 v) : p(p), v(v) {}

   /**
    * @brief Calculates distance from a point to the line
    * @param point Point to calculate distance to
    * @return Distance from the point to the line
    * @todo Implementation pending
    */
   scalar distance(const Point2 &point) {
      // TODO
   }

   /**
    * @brief Checks if this line segment intersects with another
    * @param other Other line segment to check intersection with
    * @return true if lines intersect, false otherwise
    * @note Uses a system of equations to find intersection point
    */
   bool intersect(Line2 other) {
      EqSystem2 eqsys(Mat2(v, -other.v), other.p - p);
      EqSol2 sol = eqsys.solve();
      if (sol.type == EqType::kIndependent && sol.values.x() >= 0 && sol.values.x() <= 1 &&
          sol.values.y() >= 0 && sol.values.y() <= 1)
         return true;
      else return false;
   }
};

/**
 * @brief Represents a 3D line or line segment
 */
struct Line {
   Point3 p;  ///< Point on the line
   Vec3 v;    ///< Direction vector of the line

   /**
    * @brief Constructs a line from a point and direction
    * @param p Point on the line
    * @param v Direction vector
    */
   Line(Point3 p, Vec3 v) : p(p), v(v) {}

   /**
    * @brief Calculates distance from a point to the line
    * @param point Point to calculate distance to
    * @return Distance from the point to the line
    * @todo Implementation pending
    */
   scalar distance(const Point3 &point) {
      // TODO
   }
};

/**
 * @brief Represents a plane in 3D space using the general form Ax + By + Cz + D = 0
 */
struct Plane {
   scalar A, B, C, D;  ///< Plane coefficients in Ax + By + Cz + D = 0 form

   /**
    * @brief Constructs a plane from three points
    * @param pt0 First point
    * @param pt1 Second point
    * @param pt2 Third point
    * @param ccw_winding If true, normal direction follows counter-clockwise winding rule
    */
   Plane(Point3 pt0, Point3 pt1, Point3 pt2, bool ccw_winding = false) {
      Vec3 normal = (pt2 - pt0).cross(pt1 - pt0).normalized();
      normal *= (1 - 2 * ccw_winding);  // Reverse normal if ccw winding

      A = normal.x();
      B = normal.y();
      C = normal.z();
      D = -(A * pt0.x() + B * pt0.y() + C * pt0.z());
   }

   /**
    * @brief Constructs a plane from normal vector and distance from origin
    * @param normal Normal vector to the plane
    * @param distance Signed distance from origin to plane along normal
    */
   Plane(Vec3 normal, scalar distance) : A(normal.x()), B(normal.y()), C(normal.z()), D(distance) {}

   /**
    * @brief Gets the normal vector of the plane
    * @return Normal vector (A, B, C)
    */
   Vec3 normal() const { return {A, B, C}; }

   /**
    * @brief Calculates distance from a point to the plane
    * @param pt Point to calculate distance to
    * @return Absolute distance from point to plane
    */
   scalar distance(Point3 pt) {
      return fabs(A * pt.x() + B * pt.y() + C * pt.z() + D) / sqrtf(A * A + B * B + C * C);
   }

   /**
    * @brief Determines which side of the plane a point lies on
    * @param pt Point to test
    * @return 1 if point is on positive side, -1 if on negative side
    */
   int side(Point3 pt) { return A * pt.x() + B * pt.y() + C * pt.z() - D > 0 ? 1 : -1; }

   /**
    * @brief Calculates intersection point with a line
    * @param s Line to intersect with
    * @return Point of intersection
    * @todo Implementation pending
    */
   Point3 intersection(Line s) {
      // TODO
   }

   /**
    * @brief Prints the plane equation to stdout
    */
   void print() { printf("%+.4fx %+.4fy %+.4fz %+.4f = 0\n", A, B, C, D); }
};

/**
 * @brief Represents a regular polygon or star shape in 2D
 */
struct RegularPolygon {
   Vec2 center;        ///< Center point of the polygon
   int num_vertex;     ///< Number of vertices
   scalar scale;       ///< Scale factor for size
   scalar rotation;    ///< Rotation angle in radians
   bool star;          ///< If true, creates a star shape
   scalar star_ratio;  ///< For star shapes, ratio of inner to outer radius

   /**
    * @brief Default constructor - initializes all values to zero
    */
   RegularPolygon()
       : center({0, 0}), num_vertex(0), scale(0), rotation(0), star(0), star_ratio(0) {}

   /**
    * @brief Constructs a regular polygon or star
    * @param center Center point
    * @param num_vertex Number of vertices
    * @param scale Scale factor
    * @param rotation Rotation angle in radians
    * @param star If true, creates a star shape
    * @param star_ratio For star shapes, ratio of inner to outer radius
    */
   RegularPolygon(Vec2 center, int num_vertex, scalar scale, scalar rotation = 0, bool star = false,
                  scalar star_ratio = 0.5f)
       : center(center),
         num_vertex(num_vertex),
         scale(scale),
         rotation(rotation),
         star(star),
         star_ratio(star_ratio) {};

   /**
    * @brief Generates path points for the polygon/star
    * @return Vector of point coordinates in [x1,y1,x2,y2,...] format
    */
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

   /**
    * @brief Prints the polygon parameters to stdout
    */
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

#endif /* __AMATH_GEOMETRY_H__ */