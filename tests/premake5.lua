workspace "AMathLib"
   configurations { "Debug", "Release" }
   language "C++"
   platforms { "x64" }
   architecture "x86_64"

   defines {
      "_CRT_SECURE_NO_WARNINGS",
    } 

   filter {"configurations:Debug"}
      defines { "DEBUG, _DEBUG" }
      symbols "On"
      targetsuffix "_d"
      runtime "Debug"
      symbols "Full"

   filter {"configurations:Release"}
      defines { "NDEBUG" }
      optimize "On"
      runtime "Release"


project "Tests"
   kind "ConsoleApp"
   targetdir ("bin/")
   objdir ("bin/obj/")

   files {
      "test_vec_mat.cc",
      "test_eqsolver.cc",
      "test_main.cc",
   }
   
   includedirs {
     "doctest/",
     "../include/",
   }

   


      