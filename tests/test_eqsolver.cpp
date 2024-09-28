/*
Test equation solver in 2D and 3D
Generated equations using this website:
https://courses.lumenlearning.com/wmopen-collegealgebra/chapter/introduction-systems-of-linear-equations-two-variables/
*/
#include <amath_eq.h>
#include <cassert>
#include <climits>
#include <cstdio>

using namespace amath;

int main() {
   printf("\n----------------------");
   printf("\nTESTS FOR 2D EQUATIONS");
   printf("\n----------------------");

   EqSystem2 system2;
   EqSol2 sol2;

   printf("\n\nTEST 1 - INDEPENDENT SYSTEM 1\n");
   system2 = EqSystem2(Eq2(2, -1, 2), Eq2(3, 6, 3));
   sol2 = system2.solve();

   system2.print();
   assert(sol2.type == EqType::kIndependent);
   assert(sol2.values.x() == 1 && sol2.values.y() == 0);

   printf("\n\nTEST 2 - INDEPENDENT SYSTEM 2\n");
   system2 = EqSystem2(Eq2(3, 4, -2), Eq2(-3, -1, 5));
   sol2 = system2.solve();

   system2.print();
   assert(sol2.type == EqType::kIndependent);
   assert(sol2.values.x() == -2 && sol2.values.y() == +1);

   printf("\n\nTEST 3 - DEPENDENT SYSTEM 1\n");
   system2 = EqSystem2(Eq2(3, 2, 6), Eq2(6, 4, 12));
   sol2 = system2.solve(false);

   system2.print();
   assert(sol2.type == EqType::kDependent);
   assert(sol2.values.x() == 0 && sol2.values.y() == 3);

   printf("\n\nTEST 4 - DEPENDENT SYSTEM 2\n");
   system2 = EqSystem2(Eq2(2, -4, 6), Eq2(3, -6, 9));
   sol2 = system2.solve(false);

   system2.print();
   assert(sol2.type == EqType::kDependent);
   assert(sol2.values.x() == 0 && sol2.values.y() == -1.5);

   printf("\n\nTEST 5 - INCONSISTENT SYSTEM 1\n");
   system2 = EqSystem2(Eq2(6, 3, 9), Eq2(4, 2, 12));
   sol2 = system2.solve(false);

   system2.print();
   assert(sol2.type == EqType::kInconsistent);

   printf("\n\nTEST 6 - INCONSISTENT SYSTEM 2\n");
   system2 = EqSystem2(Eq2(0, 2, 2), Eq2(0, 4, 2));
   sol2 = system2.solve(false);

   system2.print();
   assert(sol2.type == EqType::kInconsistent);

   // Test 3D equations
   printf("\n----------------------");
   printf("\nTESTS FOR 3D EQUATIONS");
   printf("\n----------------------");
   EqSystem3 system3;
   EqSol3 sol3;

   printf("\n\nTEST 1 - INDEPENDENT SYSTEM 1\n");
   system3 = EqSystem3(Eq3(1, -2, 3, 7), Eq3(-3, 1, 2, -5), Eq3(2, 2, 1, 3));
   sol3 = system3.solve();

   system3.print();
   assert(sol3.type == EqType::kIndependent);
   assert(sol3.values.x() == 2 && sol3.values.y() == -1 && sol3.values.z() == 1);

   printf("\n\nTEST 2 - INDEPENDENT SYSTEM 2\n");
   system3 = EqSystem3(Eq3(3, -2, 1, -5), Eq3(1, 3, -1, 12), Eq3(1, 1, 2, 0));
   sol3 = system3.solve();

   system3.print();
   assert(sol3.type == EqType::kIndependent);
   assert(sol3.values.x() == 1 && sol3.values.y() == 3 && sol3.values.z() == -2);

   printf("\n\nTEST 3 - INDEPENDENT SYSTEM 3\n");
   system3 = EqSystem3(Eq3(1, -1, 1, 8), Eq3(3, 3, -9, -6), Eq3(7, -2, 5, 2));
   sol3 = system3.solve(false);

   system3.print();
   assert(sol3.type == EqType::kIndependent);
   assert(sol3.values.x() == -0.625 && sol3.values.y() == -12.25 && sol3.values.z() == -3.625);

   printf("\n\nTEST 4 - DEPENDENT SYSTEM 1\n");
   system3 = EqSystem3(Eq3(1, -3, -4, 3), Eq3(3, 4, -1, 13), Eq3(2, -19, -19, 2));
   sol3 = system3.solve(false);

   system3.print();
   assert(sol3.type == EqType::kDependent);

   printf("\n\nTEST 5 - INCONSISTENT SYSTEM 1\n");
   system3 = EqSystem3(Eq3(2, -1, 1, -1), Eq3(4, 3, 5, 1), Eq3(0, 5, 3, 4));
   sol3 = system3.solve(false);

   system3.print();
   assert(sol3.type == EqType::kInconsistent);

   printf("\n\nTEST 6 - INCONSISTENT SYSTEM 2\n");
   system3 = EqSystem3(Eq3(1, 1, 1, 4), Eq3(2, -4, -1, -1), Eq3(1, -1, 0, 2));
   sol3 = system3.solve(false);

   system3.print();
   assert(sol3.type == EqType::kInconsistent);

   // Test 3D equations
   printf("\n----------------------");
   printf("\nTESTS FOR 4D EQUATIONS");
   printf("\n----------------------");
   EqSystem4 system4;
   EqSol4 sol4;

   printf("\n\nTEST 1 - INDEPENDENT SYSTEM 1\n");

   system4 = EqSystem4(
       Vec4(1, 0, -1, 0), -2, Vec4(0, 2, 0, -1), 0, Vec4(1, -2, 1, 0), 0, Vec4(0, 0, -1, 1), 1);
   sol4 = system4.solve();
   system4.print();
   assert(sol4.type == EqType::kIndependent);
   assert(sol4.values.x() == 1 && sol4.values.y() == 2 && sol4.values.z() == 3 &&
          sol4.values.w() == 4);

   printf("\n\nTEST 2 - INCONSISTENT SYSTEM 1\n");

   system4 = EqSystem4(
       Vec4(1, -1, -5, 3), -1, Vec4(1, 1, 1, -3), 0, Vec4(0, 1, 5, -3), 1, Vec4(1, -2, -10, 6), -1);
   sol4 = system4.solve(false);
   system4.print();
   assert(sol4.type == EqType::kInconsistent);
}