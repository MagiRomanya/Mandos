#ifndef PHYSICS_RENDER_H_
#define PHYSICS_RENDER_H_

#include <raylib.h>
#include <raymath.h>
#include "linear_algebra.hpp"

#include "physics_state.hpp"
#include "simulable_generator.hpp"

class MassSpringRenderer {
    public:
        MassSpringRenderer(Mesh mesh, SimulableBounds bounds)
            :mesh(mesh), dof_index(bounds.dof_index), nDoF(bounds.nDoF)
        {
            texture = LoadTexture("resources/textures/bricks.png");
        }

        ~MassSpringRenderer() {
            UnloadTexture(texture);
        }

        void draw(const PhysicsState& state) {
            Vec mesh_vertices = state.x.segment(dof_index, nDoF);
            UpdateMeshBuffer(mesh, 0, mesh_vertices.data(), mesh_vertices.size()*sizeof(Scalar), 0);
            Material material = LoadMaterialDefault();
            SetMaterialTexture(&material, MATERIAL_MAP_DIFFUSE, texture);
            DrawMesh(mesh, material, MatrixIdentity());
        }

    private:
        Texture2D texture;
        const Mesh mesh;
        const unsigned int dof_index, nDoF;
};

struct PhysicsRenderers {
    std::vector<MassSpringRenderer> mass_spring;

    void draw(const PhysicsState& state) {
        for (unsigned int i = 0; i < mass_spring.size(); i++) {
            mass_spring[i].draw(state);
        }
    }
};

#endif // PHYSICS_RENDER_H_
