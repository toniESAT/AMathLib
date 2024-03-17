#pragma once

#include "amath_core.h"

#include <fmt/format.h>
#include <vector>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265359

using amath::Vec2;

struct Polyshape {
   Vec2 center;
   int num_vertex;
   float scale;
   float rotation;
   bool star;
   float star_ratio;

   Polyshape() : center({0, 0}), num_vertex(0), scale(0), rotation(0), star(0), star_ratio(0) {}

   Polyshape(Vec2 center, int num_vertex, float scale, float rotation = 0, bool star = false,
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