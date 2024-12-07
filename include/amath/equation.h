#ifndef __AMATH_EQUATION_H__
#define __AMATH_EQUATION_H__

#include <math.h>
#include "linalg.h"
#include "utils.h"

namespace amath {

/** @brief Type alias for representing a 2D equation (ax + by = c) using Vec3 */
typedef Vec3 Eq2;
/** @brief Type alias for representing a 3D equation (ax + by + cz = d) using Vec4 */
typedef Vec4 Eq3;

/**
 * @brief Enumeration for different types of equation system solutions
 */
enum class EqType {
   kUnsolved,       ///< System hasn't been solved yet
   kIndependent,    ///< System has a unique solution
   kDependent,      ///< System has infinite solutions
   kInconsistent,   ///< System has no solutions
   kNotIndependent  ///< System is either dependent or inconsistent
};

/**
 * @brief Solution structure for 2D equation systems
 */
struct EqSol2 {
   Vec2 values;  ///< Solution values (x, y)
   EqType type;  ///< Type of solution
};

/**
 * @brief 2D linear equation system solver
 *
 * Solves systems of the form:
 * a1x + b1y = c1
 * a2x + b2y = c2
 *
 * Using Cramer's rule for independent systems.
 */
struct EqSystem2 {
   Mat2 coef;        ///< Coefficient matrix [a1 b1; a2 b2]
   Vec2 constants;   ///< Constants vector [c1; c2]
   EqSol2 solution;  ///< Solution of the system

   /**
    * @brief Constructs system from coefficient matrix and constants
    * @param coefficients Matrix of coefficients
    * @param constants Vector of constants
    */
   EqSystem2(Mat2 coefficients, Vec2 constants) : coef(coefficients), constants(constants) {};

   /**
    * @brief Default constructor - initializes all values to 0
    */
   EqSystem2() : EqSystem2(Mat2{}, Vec2{}) {};

   /**
    * @brief Constructs system from two equations in Vec3 form
    * @param eq1 First equation [a1, b1, c1]
    * @param eq2 Second equation [a2, b2, c2]
    */
   EqSystem2(const Eq2 &eq1, const Eq2 &eq2) {
      coef[0] = eq1[0];       // a1
      coef[2] = eq1[1];       // b1
      coef[1] = eq2[0];       // a2
      coef[3] = eq2[1];       // b2
      constants[0] = eq1[2];  // c1
      constants[1] = eq2[2];  // c2
   };

   /**
    * @brief Prints the equation system and its solution to stdout
    */
   void print() {
      printf("2D linear equation system:\n");

      for (int i = 0; i < 2; i++)
         printf("\tEq %d: %+.3f * x %+.3f * y = %+.3f\n",
                i + 1,
                coef[i + 2 * 0],
                coef[i + 2 * 1],
                constants[i]);

      printf("Solution: ");
      switch (solution.type) {
         case EqType::kIndependent:
            printf("x = %+.3f, y = %+.3f:\n", solution.values.x(), solution.values.y());
            break;
         case EqType::kNotIndependent:
            printf("System is not independent, no solution was calculated.");
            break;
         case EqType::kDependent:
            printf("Dependent system, an example solution is: x = %+.3f, y = %+.3f:\n",
                   solution.values.x(),
                   solution.values.y());
            break;
         case EqType::kInconsistent: printf("Inconsistent system, has no solutions\n"); break;
         case EqType::kUnsolved: printf("Unsolved\n"); break;
      }
   }

   /**
    * @brief Solves the equation system
    *
    * Uses Cramer's rule to solve the system when it's independent.
    * For dependent systems, provides an example solution.
    * For inconsistent systems, indicates no solution exists.
    *
    * @param only_independent If true, only solves for independent systems
    * @return Reference to the solution structure
    * @note Returns by reference to avoid copying
    */
   EqSol2 &solve(bool only_independent = true) {
      scalar D = coef.det();
      bool zero_determinant = almostZero(D);

      if (zero_determinant && only_independent) {
         solution.values = Vec2::nan();
         solution.type = EqType::kNotIndependent;
         return solution;
      }

      scalar D_x = coef[3] * constants[0] - coef[2] * constants[1];  // b2*c1-b1*c2
      scalar D_y = coef[0] * constants[1] - coef[1] * constants[0];  // a1*c2-a2*c1

      if (!zero_determinant) {
         solution.values = {D_x / D, D_y / D};
         solution.type = EqType::kIndependent;
      } else {
         if (almostZero(D_x) && almostZero(D_y)) {
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

/**
 * @brief Solution structure for 3D equation systems
 */
struct EqSol3 {
   Vec3 values;  ///< Solution values (x, y, z)
   EqType type;  ///< Type of solution
};

/**
 * @brief 3D linear equation system solver
 *
 * Solves systems of the form:
 * a1x + b1y + c1z = d1
 * a2x + b2y + c2z = d2
 * a3x + b3y + c3z = d3
 *
 * Using Cramer's rule for independent systems.
 */
struct EqSystem3 {
   Mat3 coef;        ///< Coefficient matrix [a1 b1 c1; a2 b2 c2; a3 b3 c3]
   Vec3 constants;   ///< Constants vector [d1; d2; d3]
   EqSol3 solution;  ///< Solution of the system

   /**
    * @brief Constructs system from coefficient matrix and constants
    * @param coefficients Matrix of coefficients
    * @param constants Vector of constants
    */
   EqSystem3(Mat3 coefficients, Vec3 constants) : coef(coefficients), constants(constants) {};

   /**
    * @brief Default constructor - initializes all values to 0
    */
   EqSystem3() : EqSystem3(Mat3{}, Vec3{}) {};

   /**
    * @brief Constructs system from three equations in Vec4 form
    * @param eq1 First equation [a1, b1, c1, d1]
    * @param eq2 Second equation [a2, b2, c2, d2]
    * @param eq3 Third equation [a3, b3, c3, d3]
    */
   EqSystem3(const Eq3 &eq1, const Eq3 &eq2, const Eq3 &eq3) {
      coef.setRow(0, {eq1.x(), eq1.y(), eq1.z()});
      coef.setRow(1, {eq2.x(), eq2.y(), eq2.z()});
      coef.setRow(2, {eq3.x(), eq3.y(), eq3.z()});

      constants[0] = eq1.w();
      constants[1] = eq2.w();
      constants[2] = eq3.w();
   };

   /**
    * @brief Prints the equation system and its solution to stdout
    */
   void print() {
      printf("3D linear equation system:\n");

      for (int i = 0; i < 3; i++)
         printf("\tEq %d: %+.3f * x %+.3f * y %+.3f * z = %+.3f\n",
                i + 1,
                coef[i + 3 * 0],
                coef[i + 3 * 1],
                coef[i + 3 * 2],
                constants[i]);

      printf("Solution: ");
      switch (solution.type) {
         case EqType::kIndependent:
            printf("x = %+.3f, y = %+.3f, z = %+.3f\n",
                   solution.values.x(),
                   solution.values.y(),
                   solution.values.z());
            break;
         case EqType::kNotIndependent:
            printf("System is not independent, no solution was calculated.");
            break;
         case EqType::kDependent: printf("Dependent system, infinite solutions\n"); break;
         case EqType::kInconsistent: printf("Inconsistent system, has no solutions\n"); break;
         case EqType::kUnsolved: printf("Unsolved\n"); break;
      }
   }

   /**
    * @brief Solves the equation system
    *
    * Uses Cramer's rule to solve the system when it's independent.
    * For dependent and inconsistent systems, indicates the solution type
    * but does not compute example solutions.
    *
    * @param only_independent If true, only solves for independent systems
    * @return Reference to the solution structure
    * @note Returns by reference to avoid copying
    */
   EqSol3 &solve(bool only_independent = true) {
      scalar D = coef.det();
      bool zero_determinant = almostZero(D);

      if (zero_determinant && only_independent) {
         solution.values = Vec3::nan();
         solution.type = EqType::kNotIndependent;
         return solution;
      }

      // Calculate determinants for Cramer's rule
      scalar D_x = constants[0] * coef[4] * coef[8] + constants[2] * coef[3] * coef[7] +
                   constants[1] * coef[5] * coef[6] - constants[2] * coef[4] * coef[6] -
                   constants[1] * coef[3] * coef[8] - constants[0] * coef[5] * coef[7];

      scalar D_y = coef[0] * constants[1] * coef[8] + coef[2] * constants[0] * coef[7] +
                   coef[1] * constants[2] * coef[6] - coef[2] * constants[1] * coef[6] -
                   coef[1] * constants[0] * coef[8] - coef[0] * constants[2] * coef[7];

      scalar D_z = coef[0] * coef[4] * constants[2] + coef[2] * coef[3] * constants[1] +
                   coef[1] * coef[5] * constants[0] - coef[2] * coef[4] * constants[0] -
                   coef[1] * coef[3] * constants[2] - coef[0] * coef[5] * constants[1];

      if (!zero_determinant) {
         solution.values = {D_x / D, D_y / D, D_z / D};
         solution.type = EqType::kIndependent;
      } else {
         solution.values = Vec3::nan();
         if (almostZero(D_x) && almostZero(D_y) && almostZero(D_z)) {
            solution.type = EqType::kDependent;
         } else {
            solution.type = EqType::kInconsistent;
         }
      }

      return solution;
   }
};

/**
 * @brief Solution structure for 4D equation systems
 */
struct EqSol4 {
   Vec4 values;  ///< Solution values (x, y, z, w)
   EqType type;  ///< Type of solution
};

/**
 * @brief 4D linear equation system solver
 *
 * Solves systems of the form:
 * a1x + b1y + c1z + d1w = k1
 * a2x + b2y + c2z + d2w = k2
 * a3x + b3y + c3z + d3w = k3
 * a4x + b4y + c4z + d4w = k4
 *
 * Using Cramer's rule for independent systems.
 */
struct EqSystem4 {
   Mat4 coef;        ///< Coefficient matrix [a1 b1 c1 d1; a2 b2 c2 d2; a3 b3 c3 d3; a4 b4 c4 d4]
   Vec4 constants;   ///< Constants vector [k1; k2; k3; k4]
   EqSol4 solution;  ///< Solution of the system

   /**
    * @brief Constructs system from coefficient matrix and constants
    * @param coefficients Matrix of coefficients
    * @param constants Vector of constants
    */
   EqSystem4(Mat4 coefficients, Vec4 constants) : coef(coefficients), constants(constants) {};

   /**
    * @brief Default constructor - initializes all values to 0
    */
   EqSystem4() : EqSystem4(Mat4{}, Vec4{}) {};

   /**
    * @brief Constructs system from four coefficient vectors and constants
    * @param coef1 First equation coefficients [a1, b1, c1, d1]
    * @param const1 First equation constant k1
    * @param coef2 Second equation coefficients [a2, b2, c2, d2]
    * @param const2 Second equation constant k2
    * @param coef3 Third equation coefficients [a3, b3, c3, d3]
    * @param const3 Third equation constant k3
    * @param coef4 Fourth equation coefficients [a4, b4, c4, d4]
    * @param const4 Fourth equation constant k4
    */
   EqSystem4(const Vec4 &coef1, const scalar const1, const Vec4 &coef2, const scalar const2,
             const Vec4 &coef3, const scalar const3, const Vec4 &coef4, const scalar const4) {
      coef.setRow(0, coef1);
      coef.setRow(1, coef2);
      coef.setRow(2, coef3);
      coef.setRow(3, coef4);

      constants[0] = const1;
      constants[1] = const2;
      constants[2] = const3;
      constants[3] = const4;
   };

   /**
    * @brief Prints the equation system and its solution to stdout
    */
   void print() {
      printf("4D linear equation system:\n");

      for (int i = 0; i < 4; i++)
         printf("\tEq %d: %+.3f * x %+.3f * y %+.3f * z %+.3f * w = %+.3f\n",
                i + 1,
                coef[i + 4 * 0],
                coef[i + 4 * 1],
                coef[i + 4 * 2],
                coef[i + 4 * 3],
                constants[i]);

      printf("Solution: ");
      switch (solution.type) {
         case EqType::kIndependent:
            printf("x = %+.3f, y = %+.3f, z = %+.3f, w = %+.3f\n",
                   solution.values.x(),
                   solution.values.y(),
                   solution.values.z(),
                   solution.values.w());
            break;
         case EqType::kNotIndependent:
            printf("System is not independent, no solution was calculated.");
            break;
         case EqType::kDependent: printf("Dependent system, infinite solutions\n"); break;
         case EqType::kInconsistent: printf("Inconsistent system, has no solutions\n"); break;
         case EqType::kUnsolved: printf("Unsolved\n"); break;
      }
   }

   /**
    * @brief Solves the equation system
    *
    * Uses Cramer's rule to solve the system when it's independent.
    * For dependent and inconsistent systems, indicates the solution type
    * but does not compute example solutions.
    *
    * @param only_independent If true, only solves for independent systems
    * @return Reference to the solution structure
    * @note Returns by reference to avoid copying
    */
   EqSol4 &solve(bool only_independent = true) {
      scalar D = coef.det();
      bool zero_determinant = almostZero(D);

      if (zero_determinant && only_independent) {
         solution.values = Vec4::nan();
         solution.type = EqType::kNotIndependent;
         return solution;
      }

      // Use Cramer's rule with temporary matrices
      Mat4 cramer_matrix;
      Vec4 cramer_determinants;

      // Calculate determinant for each variable
      for (int i = 0; i < 4; i++) {
         cramer_matrix = coef;
         cramer_matrix.setCol(i, constants);
         cramer_determinants[i] = cramer_matrix.det();
      }

      if (!zero_determinant) {
         solution.values = cramer_determinants / D;
         solution.type = EqType::kIndependent;
      } else {
         bool is_dependent = true;
         for (int i = 0; i < 4; i++)
            is_dependent = is_dependent && almostZero(cramer_determinants[i]);

         solution.values = Vec4::nan();
         if (is_dependent) {
            solution.type = EqType::kDependent;
         } else {
            solution.type = EqType::kInconsistent;
         }
      }

      return solution;
   }
};

}  // namespace amath

#endif /* __AMATH_EQUATION_H__ */
