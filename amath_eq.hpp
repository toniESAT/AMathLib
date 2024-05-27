#pragma once

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <math.h>
#include <amath_core.hpp>
#include <amath_utils.hpp>

namespace amath {

typedef Vec3 Eq2; // For readability

enum class EqType {
   kUnsolved,
   kIndependent, // Unique solution
   kDependent,
   kInconsistent,
   kNotIndependent // Either dependent or inconsistent
};

struct EqSol2 {
   Vec2 values;
   EqType type;
};

struct EqSystem2 {
   Mat2 coef;
   Vec2 constants;
   EqSol2 solution = {Vec2::nan(), EqType::kUnsolved};

   EqSystem2(Mat2 coefficients, Vec2 constants) : coef(coefficients), constants(constants){};
   EqSystem2() : EqSystem2(Mat2{}, Vec2{}){}; // Init all to 0
   EqSystem2(const Eq2 &eq1, const Eq2 &eq2) {
      /* Current implementation:
                      Mat2 coef     Vec2 consts
         Vec3 eq1 [a1 * x + b1 * y ] = [c1]
         Vec3 eq2 [a2 * x + b2 * y ] = [c2]
      */
      coef[0] = eq1[0];
      coef[2] = eq1[1];

      coef[1] = eq2[0];
      coef[3] = eq2[1];

      constants[0] = eq1[2];
      constants[1] = eq2[2];
   };

   void print() {
      fmt::print("2D linear equation system:\n");

      for (int i = 0; i < 2; i++)
         fmt::print("\tEq {:d}: {:+.3f} * x {:+.3f} * y = {:+.3f}\n",
                    i + 1,
                    coef[i + 2 * 0],
                    coef[i + 2 * 1],
                    constants[i]);

      fmt::print("Solution: ");
      switch (solution.type) {
      case EqType::kIndependent:
         fmt::print("x = {:+.3f}, y = {:+.3f}:\n", solution.values.x(), solution.values.y());
         break;
      case EqType::kNotIndependent:
         fmt::print("System is not independent, no solution was calculated.");
         break;
      case EqType::kDependent:
         fmt::print("Dependent system, an example solution is: x = {:+.3f}, y = {:+.3f}:\n",
                    solution.values.x(),
                    solution.values.y());
         break;
      case EqType::kInconsistent: fmt::print("Inconsistent system, has no solutions\n"); break;
      case EqType::kUnsolved: fmt::print("Unsolved\n"); break;
      }
   }

   // CAUTION: return by reference
   EqSol2 &solve(bool only_independent = true) {

      float D = coef.determinant();
      bool zero_determinant = is_almost_zero(D);

      /* When  D=0,
      D_x==0   and   D_y==0 -> Dependent System
      D_x!=0   or    D_y!=0 -> Inconsistent System */

      if (zero_determinant && only_independent) {
         solution.values = Vec2::nan();
         solution.type = EqType::kNotIndependent;
         return solution; // If no unique solution, return early
      }

      float D_x = coef[3] * constants[0] - coef[2] * constants[1]; // b2*c1-b1*c2
      float D_y = coef[0] * constants[1] - coef[1] * constants[0]; // a1*c2-a2*c1

      if (!zero_determinant) {
         solution.values = {D_x / D, D_y / D};
         solution.type = EqType::kIndependent;
      } else {
         if (is_almost_zero(D_x) && is_almost_zero(D_y)) {
            // Return solution for x = 0, for example
            solution.values = {0, constants[0] / coef[2]};
            solution.type = EqType::kDependent;
         } else {
            solution.values = Vec2::nan();
            solution.type = EqType::kInconsistent;
         }
      }

      return solution;
   }
};

typedef Vec4 Eq3;

struct EqSol3 {
   Vec3 values;
   EqType type;
};

struct EqSystem3 {
   Mat3 coef;
   Vec3 constants;
   EqSol3 solution = {Vec3::nan(), EqType::kUnsolved};

   EqSystem3(Mat3 coefficients, Vec3 constants) : coef(coefficients), constants(constants){};
   EqSystem3() : EqSystem3(Mat3{}, Vec3{}){}; // Init all to 0
   EqSystem3(const Eq3 &eq1, const Eq3 &eq2, const Eq3 &eq3) {
      /* Current implementation:
                               Mat3 coef         Vec3 consts
         Vec4/Eq3 eq1 [a1 * x + b1 * y + c1 * z] = [d1]
         Vec4/Eq3 eq2 [a2 * x + b2 * y + c2 * z] = [d2]
         Vec4/Eq3 eq3 [a3 * x + b3 * y + c3 * z] = [d3]
      */

      coef.set_row(0, {eq1.x(), eq1.y(), eq1.z()});
      coef.set_row(1, {eq2.x(), eq2.y(), eq2.z()});
      coef.set_row(2, {eq3.x(), eq3.y(), eq3.z()});

      constants[0] = eq1.w();
      constants[1] = eq2.w();
      constants[2] = eq3.w();
   };

   void print() {
      fmt::print("3D linear equation system:\n");

      for (int i = 0; i < 3; i++)
         fmt::print("\tEq {:d}: {:+.3f} * x {:+.3f} * y {:+.3f} * z = {:+.3f}\n",
                    i + 1,
                    coef[i + 3 * 0],
                    coef[i + 3 * 1],
                    coef[i + 3 * 2],
                    constants[i]);
      fmt::print("Solution: ");
      switch (solution.type) {
      case EqType::kIndependent:
         fmt::print("x = {:+.3f}, y = {:+.3f}, z = {:+.3f}\n",
                    solution.values.x(),
                    solution.values.y(),
                    solution.values.z());
         break;
      case EqType::kNotIndependent:
         fmt::print("System is not independent, no solution was calculated.");
         break;
      case EqType::kDependent: fmt::print("Dependent system, infinite solutions\n"); break;
      case EqType::kInconsistent: fmt::print("Inconsistent system, has no solutions\n"); break;
      case EqType::kUnsolved: fmt::print("Unsolved\n"); break;
      }
   }
   // CAUTION: return by reference
   EqSol3 &solve(bool only_independent = true) {

      float D = coef.determinant();
      bool zero_determinant = is_almost_zero(D);

      /* When  D=0,
      D_x==0   and   D_y==0 and D_z==0-> Dependent System
      D_x!=0   or    D_y!=0 or    D_z!=0-> Inconsistent System */

      if (zero_determinant && only_independent) {
         solution.values = Vec3::nan();
         solution.type = EqType::kNotIndependent;
         return solution; // If no unique solution, return early
      }

      float D_x = constants[0] * coef[4] * coef[8] + constants[2] * coef[3] * coef[7] +
                  constants[1] * coef[5] * coef[6] - constants[2] * coef[4] * coef[6] -
                  constants[1] * coef[3] * coef[8] - constants[0] * coef[5] * coef[7];

      float D_y = coef[0] * constants[1] * coef[8] + coef[2] * constants[0] * coef[7] +
                  coef[1] * constants[2] * coef[6] - coef[2] * constants[1] * coef[6] -
                  coef[1] * constants[0] * coef[8] - coef[0] * constants[2] * coef[7];

      float D_z = coef[0] * coef[4] * constants[2] + coef[2] * coef[3] * constants[1] +
                  coef[1] * coef[5] * constants[0] - coef[2] * coef[4] * constants[0] -
                  coef[1] * coef[3] * constants[2] - coef[0] * coef[5] * constants[1];

      if (!zero_determinant) {
         solution.values = {D_x / D, D_y / D, D_z / D};
         solution.type = EqType::kIndependent;
      } else {
         solution.values = Vec3::nan();
         if (is_almost_zero(D_x) && is_almost_zero(D_y) && is_almost_zero(D_z)) {
            // Could be dependent in several different manners
            // Currently, don't bother calculating a possible solution, just return NaN
            // solution.values = Vec3::nan();
            solution.type = EqType::kDependent;
         } else {
            // solution.values = Vec3::nan();
            solution.type = EqType::kInconsistent;
         }
      }

      return solution;
   }
};

struct EqSol4 {
   Vec4 values;
   EqType type;
};

struct EqSystem4 {
   Mat4 coef;
   Vec4 constants;
   EqSol4 solution = {Vec4::nan(), EqType::kUnsolved};

   EqSystem4(Mat4 coefficients, Vec4 constants) : coef(coefficients), constants(constants){};
   EqSystem4() : EqSystem4(Mat4{}, Vec4{}){}; // Init all to 0
   EqSystem4(const Vec4 &coef1, const float const1, const Vec4 &coef2, const float const2,
             const Vec4 &coef3, const float const3, const Vec4 &coef4, const float const4) {

      /*      args           Mat4 coef                   args   Vec4 consts
         Vec4 coef1 [a1 * x + b1 * y + c1 * z + d1 * w] =  (float const1)   [k1]
         Vec4 coef2 [a2 * x + b2 * y + c2 * z + d2 * w] =  (float const2)   [k2]
         Vec4 coef3 [a3 * x + b3 * y + c3 * z + d3 * w] =  (float const3)   [k3]
         Vec4 coef3 [a3 * x + b3 * y + c3 * z + d3 * w] =  (float const4)   [k3]
      */

      coef.set_row(0, coef1);
      coef.set_row(1, coef2);
      coef.set_row(2, coef3);
      coef.set_row(3, coef4);

      constants[0] = const1;
      constants[1] = const2;
      constants[2] = const3;
      constants[3] = const4;
   };

   void print() {
      fmt::print("4D linear equation system:\n");

      for (int i = 0; i < 4; i++)
         fmt::print("\tEq {:d}: {:+.3f} * x {:+.3f} * y {:+.3f} * z {:+.3f} * w = {:+.3f}\n",
                    i + 1,
                    coef[i + 4 * 0],
                    coef[i + 4 * 1],
                    coef[i + 4 * 2],
                    coef[i + 4 * 3],
                    constants[i]);
      fmt::print("Solution: ");
      switch (solution.type) {
      case EqType::kIndependent:
         fmt::print("x = {:+.3f}, y = {:+.3f}, z = {:+.3f},  w = {:+.3f}\n",
                    solution.values.x(),
                    solution.values.y(),
                    solution.values.z(),
                    solution.values.w());
         break;
      case EqType::kNotIndependent:
         fmt::print("System is not independent, no solution was calculated.");
         break;
      case EqType::kDependent: fmt::print("Dependent system, infinite solutions\n"); break;
      case EqType::kInconsistent: fmt::print("Inconsistent system, has no solutions\n"); break;
      case EqType::kUnsolved: fmt::print("Unsolved\n"); break;
      }
   }
   // CAUTION: return by reference
   EqSol4 &solve(bool only_independent = true) {

      float D = coef.determinant();
      bool zero_determinant = is_almost_zero(D);

      /* When  D=0,
      D_x==0   and   D_y==0 and D_z==0 and D_w==0-> Dependent System
      D_x!=0   or    D_y!=0 or    D_z!=0 or  D_w!=0-> Inconsistent System */

      if (zero_determinant && only_independent) {
         solution.values = Vec4::nan();
         solution.type = EqType::kNotIndependent;
         return solution; // If no unique solution, return early
      }

      Mat4 cramer_matrix;
      Vec4 cramer_determinants;
      for (int i = 0; i < 4; i++) {
         cramer_matrix = coef;
         cramer_matrix.set_column(i, constants);
         cramer_determinants[i] = cramer_matrix.determinant();
      }

      if (!zero_determinant) {
         solution.values = cramer_determinants / D;
         solution.type = EqType::kIndependent;
      } else {
         bool is_dependent = true;
         for (int i = 0; i < 4; i++)
            is_dependent = is_dependent && is_almost_zero(cramer_determinants[i]);

         solution.values = Vec4::nan();
         if (is_dependent) {
            // Could be dependent in several different manners
            // Currently, don't bother calculating a possible solution, just return NaN
            // solution.values = Vec3::nan();
            solution.type = EqType::kDependent;
         } else {
            // solution.values = Vec3::nan();
            solution.type = EqType::kInconsistent;
         }
      }

      return solution;
   }
};

} // namespace amath