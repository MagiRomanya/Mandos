#ifndef DRAW_VECTOR_H__
#define DRAW_VECTOR_H__

#include "raylib.h"
#include "linear_algebra.hpp"
#include "utility_functions.hpp"
#include <iostream>

inline void DrawVector(Vec3 position, Vec3 vector, Color color, Scalar scale = 1) {
  const Vec3 v = (vector - position).normalized();
  const Vec3 up = Vec3(0.0, 1.0, 0.0);
  Vec3 rotation_axis = - skew(v) * up;
  if (rotation_axis == Vec3::Zero()) {
    rotation_axis = Vec3(1,0,0);
  }
  const Scalar angle = std::acos(v.y());
  const Mat3 rot = compute_rotation_matrix_rodrigues(rotation_axis.normalized()*angle);
  // std::cout << "Angle " << angle << std::endl;
  // std::cout << "Axis " << rotation_axis << std::endl;
  // std::cout << rot << std::endl;
  Matrix transform = {scale*rot(0,0), scale*rot(0,1), scale*rot(0,2), position.x(),
                      scale*rot(1,0), scale*rot(1,1), scale*rot(1,2), position.y(),
                      scale*rot(2,0), scale*rot(2,1), scale*rot(2,2), position.z(),
                           0.0,      0.0,      0.0,          1.0};
  // render
  const Mesh vector_mesh = LoadMeshTinyOBJ("img/obj/vector.obj");
  Material material = LoadMaterialDefault();
  Image im_color = GenImageColor(1, 1, color);
  Texture2D texture = LoadTextureFromImage(im_color);
  SetMaterialTexture(&material, MATERIAL_MAP_DIFFUSE, texture);
  DrawMesh(vector_mesh, material, transform);
  UnloadMesh(vector_mesh);
  UnloadTexture(texture);
}
#endif // DRAW_VECTOR_H__
