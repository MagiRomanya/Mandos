#ifndef PHYSICS_RENDER_H_
#define PHYSICS_RENDER_H_

#include "linear_algebra.hpp"
#include <raylib.h>
#include <raymath.h>

#include "physics_state.hpp"

class PhysicsRender {
    public:
        PhysicsRender(Mesh mesh, unsigned int dof_index, unsigned int nDoF, Color color)
            : dof_index(dof_index), nDoF(nDoF), color(color), mesh(mesh) {}

        virtual void draw(const PhysicsState& state) = 0;

    protected:
        const Mesh mesh;
        const unsigned int dof_index, nDoF;
        const Color color;
};

class MassSpringRender : public PhysicsRender {
    public:
        MassSpringRender(Mesh mesh, unsigned int dof_index, unsigned int nDoF, Color color)
            : PhysicsRender(mesh, dof_index, nDoF, color) {}

        void draw(const PhysicsState& state) override {
            Vec mesh_vertices = state.x.segment(dof_index, dof_index + nDoF);
            UpdateMeshBuffer(mesh, 0, mesh_vertices.data(), mesh_vertices.size()*sizeof(Scalar), 0);
            Material material = LoadMaterialDefault();
            Texture2D texture = LoadTexture("img/textures/bricks.png");
            SetMaterialTexture(&material, MATERIAL_MAP_DIFFUSE, texture);
            DrawMesh(mesh, material, MatrixIdentity());
        }
};

#endif // PHYSICS_RENDER_H_
