#include <cstring>
#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>

#include "inertia_energies.hpp"
#include "linear_algebra.hpp"
#include "mandos.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "raymath.h"
#include "render/simulation_visualization.hpp"
#include "rigid_body.hpp"
#include "simulable_generator.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"


inline Matrix matrix_eigen_to_raylib(const Mat4& m) {
    Matrix r = {
    m(0,0), m(0,1), m(0, 2), m(0, 3),
    m(1,0), m(1,1), m(1, 2), m(1, 3),
    m(2,0), m(2,1), m(2, 2), m(2, 3),
    m(3,0), m(3,1), m(3, 2), m(3, 3),
    };
    return r;
}

int main(void) {

    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos simulator");
    Camera3D camera = create_camera();
    //--------------------------------------------------------------------------------------

    Simulation simulation;
    const Scalar RB_MASS = 1;

    const Mesh RB_mesh = GenMeshPlane(4.0, 1.0, 20, 20);
    const unsigned int n_vertices = RB_mesh.vertexCount * 3;
    const unsigned int n_indices = RB_mesh.triangleCount * 3;
    const std::vector<Scalar> vertices(RB_mesh.vertices, RB_mesh.vertices + n_vertices);
    const std::vector<unsigned int> indices(RB_mesh.indices, RB_mesh.indices + n_indices);
    const Mat3 inertia_tensor = compute_initial_inertia_tensor_PARTICLES(RB_MASS, vertices);

    const RigidBodyHandle rb1 = RigidBodyHandle(simulation, RB_MASS, inertia_tensor)
        .set_COM_initial_position(simulation, Vec3(0,0,0))
        .add_gravity(simulation, -1);

    const RigidBodyHandle rb2 = RigidBodyHandle(simulation, RB_MASS, inertia_tensor)
        .set_COM_initial_position(simulation, Vec3(4,0,0))
        .add_gravity(simulation, -1);

    const RigidBodyHandle rb3 = RigidBodyHandle(simulation, RB_MASS, inertia_tensor)
        .set_COM_initial_position(simulation, Vec3(8,0,0))
        .add_gravity(simulation, -1);

    simulation.constraints.rb_point_constraints.emplace_back(0, rb1.rb, rb2.rb, Vec3(2,0,0), Vec3(-2,0,0));
    simulation.constraints.rb_point_constraints.emplace_back(3, rb2.rb, rb3.rb, Vec3(2,0,0), Vec3(-2,0,0));

    simulation.frozen_dof.push_back(0);
    simulation.frozen_dof.push_back(1);
    simulation.frozen_dof.push_back(2);

    SetTargetFPS(200);
    //--------------------------------------------------------------------------------------

    bool simulation_paused = true;
    PhysicsState state = simulation.initial_state;

    Material material1 = LoadMaterialDefault();
    material1.maps[MATERIAL_MAP_DIFFUSE].color = RED;
    Material material2 = LoadMaterialDefault();
    material2.maps[MATERIAL_MAP_DIFFUSE].color = GREEN;
    Material material3 = LoadMaterialDefault();
    material3.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;

    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera, CAMERA_ORBITAL);
        /// Keyboard controls
        if (IsKeyPressed(KEY_Q)) break;

        if (IsKeyPressed(KEY_H)) simulation_paused = !simulation_paused;
        if (IsKeyPressed(KEY_R)) state = simulation.initial_state;
        if (!simulation_paused or IsKeyPressed(KEY_G)) {
            EnergyAndDerivatives f(0);
            simulation_step(simulation, state, f);
            std::cout << "Energy " << f.energy << std::endl;
        }
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);
            {
                DrawGrid(30, 1.0f);
                // simulation_render_simulables_and_energies(simulation, state);
                Matrix rb1_m = matrix_eigen_to_raylib(rb1.get_RB_transformation(state));
                Matrix rb2_m = matrix_eigen_to_raylib(rb2.get_RB_transformation(state));
                Matrix rb3_m = matrix_eigen_to_raylib(rb3.get_RB_transformation(state));
                DrawMesh(RB_mesh, material1, rb1_m);
                DrawMesh(RB_mesh, material2, rb2_m);
                DrawMesh(RB_mesh, material3, rb3_m);
            }
            EndMode3D();

            DrawFPS(10, 10);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }
    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
