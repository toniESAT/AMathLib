/*
Test equation solver in 2D and 3D
Generated equations using this website:
https://courses.lumenlearning.com/wmopen-collegealgebra/chapter/introduction-systems-of-linear-equations-two-variables/
*/
#include "amath/equation.h"

#include "doctest.h"

TEST_CASE("TESTS FOR 2D EQUATIONS") {
   using amath::Eq2;
   using amath::EqSol2;
   using amath::EqSystem2;
   using amath::EqType;

   EqSystem2 system2;
   EqSol2 sol2;

   SUBCASE("TEST 1 - INDEPENDENT SYSTEM 1") {
      system2 = EqSystem2(Eq2(2, -1, 2), Eq2(3, 6, 3));
      sol2 = system2.solve();

      CHECK(sol2.type == amath::EqType::kIndependent);
      CHECK(sol2.values.x() == 1);
      CHECK(sol2.values.y() == 0);
   }

   SUBCASE("TEST 2 - INDEPENDENT SYSTEM 2") {
      system2 = EqSystem2(Eq2(3, 4, -2), Eq2(-3, -1, 5));
      sol2 = system2.solve();

      CHECK(sol2.type == EqType::kIndependent);
      CHECK(sol2.values.x() == -2);
      CHECK(sol2.values.y() == +1);
   }

   SUBCASE("TEST 3 - DEPENDENT SYSTEM 1") {
      system2 = EqSystem2(Eq2(3, 2, 6), Eq2(6, 4, 12));
      sol2 = system2.solve(false);

      CHECK(sol2.type == EqType::kDependent);
      CHECK(sol2.values.x() == 0);
      CHECK(sol2.values.y() == 3);
   }

   SUBCASE("TEST 4 - DEPENDENT SYSTEM 2") {
      system2 = EqSystem2(Eq2(2, -4, 6), Eq2(3, -6, 9));
      sol2 = system2.solve(false);

      CHECK(sol2.type == EqType::kDependent);
      CHECK(sol2.values.x() == 0);
      CHECK(sol2.values.y() == -1.5);
   }

   SUBCASE("TEST 5 - INCONSISTENT SYSTEM 1") {
      system2 = EqSystem2(Eq2(6, 3, 9), Eq2(4, 2, 12));
      sol2 = system2.solve(false);

      CHECK(sol2.type == EqType::kInconsistent);
   }

   SUBCASE("TEST 6 - INCONSISTENT SYSTEM 2") {
      system2 = EqSystem2(Eq2(0, 2, 2), Eq2(0, 4, 2));
      sol2 = system2.solve(false);

      CHECK(sol2.type == EqType::kInconsistent);
   }
}

// Test 3D equations
TEST_CASE("\nTESTS FOR 3D EQUATIONS") {
   using amath::Eq3;
   using amath::EqSol3;
   using amath::EqSystem3;
   using amath::EqType;

   EqSystem3 system3;
   EqSol3 sol3;

   SUBCASE("TEST 1 - INDEPENDENT SYSTEM 1") {
      system3 = EqSystem3(Eq3(1, -2, 3, 7), Eq3(-3, 1, 2, -5), Eq3(2, 2, 1, 3));
      sol3 = system3.solve();

      CHECK(sol3.type == EqType::kIndependent);
      CHECK(sol3.values.x() == 2);
      CHECK(sol3.values.y() == -1);
      CHECK(sol3.values.z() == 1);
   }

   SUBCASE("TEST 2 - INDEPENDENT SYSTEM 2") {
      system3 = EqSystem3(Eq3(3, -2, 1, -5), Eq3(1, 3, -1, 12), Eq3(1, 1, 2, 0));
      sol3 = system3.solve();

      CHECK(sol3.type == EqType::kIndependent);
      CHECK(sol3.values.x() == 1);
      CHECK(sol3.values.y() == 3);
      CHECK(sol3.values.z() == -2);
   }

   SUBCASE("TEST 3 - INDEPENDENT SYSTEM 3") {
      system3 = EqSystem3(Eq3(1, -1, 1, 8), Eq3(3, 3, -9, -6), Eq3(7, -2, 5, 2));
      sol3 = system3.solve(false);

      CHECK(sol3.type == EqType::kIndependent);
      CHECK(sol3.values.x() == -0.625);
      CHECK(sol3.values.y() == -12.25);
      CHECK(sol3.values.z() == -3.625);
   }

   SUBCASE("TEST 4 - DEPENDENT SYSTEM 1") {
      system3 = EqSystem3(Eq3(1, -3, -4, 3), Eq3(3, 4, -1, 13), Eq3(2, -19, -19, 2));
      sol3 = system3.solve(false);

      CHECK(sol3.type == EqType::kDependent);
   }
   SUBCASE("TEST 5 - INCONSISTENT SYSTEM 1") {
      system3 = EqSystem3(Eq3(2, -1, 1, -1), Eq3(4, 3, 5, 1), Eq3(0, 5, 3, 4));
      sol3 = system3.solve(false);

      CHECK(sol3.type == EqType::kInconsistent);
   }

   SUBCASE("TEST 6 - INCONSISTENT SYSTEM 2") {
      system3 = EqSystem3(Eq3(1, 1, 1, 4), Eq3(2, -4, -1, -1), Eq3(1, -1, 0, 2));
      sol3 = system3.solve(false);

      CHECK(sol3.type == EqType::kInconsistent);
   }
}

// Test 3D equations
TEST_CASE("\nTESTS FOR 4D EQUATIONS") {
   using amath::EqSol4;
   using amath::EqSystem4;
   using amath::EqType;
   using amath::Vec4;

   EqSystem4 system4;
   EqSol4 sol4;

   SUBCASE("TEST 1 - INDEPENDENT SYSTEM 1") {
      system4 = EqSystem4(
          Vec4(1, 0, -1, 0), -2, Vec4(0, 2, 0, -1), 0, Vec4(1, -2, 1, 0), 0, Vec4(0, 0, -1, 1), 1);
      sol4 = system4.solve();

      CHECK(sol4.type == EqType::kIndependent);
      CHECK(sol4.values.x() == 1);
      CHECK(sol4.values.y() == 2);
      CHECK(sol4.values.z() == 3);
      CHECK(sol4.values.w() == 4);
   }

   SUBCASE("TEST 2 - INCONSISTENT SYSTEM 1") {
      system4 = EqSystem4(Vec4(1, -1, -5, 3),
                          -1,
                          Vec4(1, 1, 1, -3),
                          0,
                          Vec4(0, 1, 5, -3),
                          1,
                          Vec4(1, -2, -10, 6),
                          -1);
      sol4 = system4.solve(false);

      CHECK(sol4.type == EqType::kInconsistent);
   }
}
