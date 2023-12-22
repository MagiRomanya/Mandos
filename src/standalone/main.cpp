#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>

#include "linear_algebra.hpp"
#include "mandos.hpp"
#include "mesh.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "physics_state.hpp"
#include "raymath.h"
#include "render/simulation_visualization.hpp"
#include "rigid_body.hpp"
#include "clock.hpp"
#include "utility_functions.hpp"


int main(void) {
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos simulator");
    Camera3D camera = create_camera();
    //--------------------------------------------------------------------------------------

    Simulation simulation;
    const Scalar MASS = 1;
    const Scalar GRAVITY = -1;

    const Mesh RB_mesh = GenMeshPlane(4.0, 1.0, 10, 10);
    const unsigned int n_vertices = RB_mesh.vertexCount * 3;
    const unsigned int n_indices = RB_mesh.triangleCount * 3;
    const std::vector<Scalar> vertices(RB_mesh.vertices, RB_mesh.vertices + n_vertices);
    const std::vector<unsigned int> indices(RB_mesh.indices, RB_mesh.indices + n_indices);
    const Mat3 inertia_tensor = compute_initial_inertia_tensor_PARTICLES(MASS, vertices);

    const RigidBodyHandle rb1 = RigidBodyHandle(simulation, MASS, inertia_tensor)
        .add_gravity(GRAVITY)
        .set_COM_initial_position(Vec3(0,0,0));

    const ParticleHandle p1 = ParticleHandle(simulation, MASS)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(1,1,0));

    const ParticleHandle p2 = ParticleHandle(simulation, MASS)
        .freeze()
        .set_initial_position(Vec3(0,0,1));

    const ParticleHandle p3 = ParticleHandle(simulation, MASS)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(0,1,-1));

    const ParticleHandle p4 = ParticleHandle(simulation, MASS)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(0,1,1));

    const ParticleHandle p5 = ParticleHandle(simulation, MASS)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(1,1,1));

    const ParticleHandle p6 = ParticleHandle(simulation, MASS)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(2,1,1));

    std::vector<ParticleRigidBodyCopuling> copulings_vec;
    ParticleRigidBodyCopuling copuling = ParticleRigidBodyCopuling(rb1.rb, p1.particle, p1.particle.get_position(simulation.initial_state.x));
    copulings_vec.push_back(copuling);
    ParticleRigidBodyCopuling copuling2 = ParticleRigidBodyCopuling(rb1.rb, p3.particle, p3.particle.get_position(simulation.initial_state.x));
    copulings_vec.push_back(copuling2);
    ParticleRigidBodyCopuling copuling3 = ParticleRigidBodyCopuling(rb1.rb, p4.particle, p4.particle.get_position(simulation.initial_state.x));
    copulings_vec.push_back(copuling3);

    simulation.copulings = Copulings(copulings_vec);

    // std::cout << rb1.rb.get_COM_position(simulation.initial_state.x).transpose() << std::endl;
    // std::cout << rb1.rb.get_COM_position(simulation.initial_state.x_old).transpose() << std::endl;
    // std::cout << simulation.initial_state.x.transpose() << std::endl;

    join_particles_with_spring(simulation, p1, p2, 10.0, 1);
    join_particles_with_spring(simulation, p4, p5, 10.0, 1);
    join_particles_with_spring(simulation, p5, p6, 10.0, 1);

    SetTargetFPS(200);
    //--------------------------------------------------------------------------------------

    bool simulation_paused = true;
    PhysicsState state = simulation.initial_state;

    Material material1 = LoadMaterialDefault();
    material1.maps[MATERIAL_MAP_DIFFUSE].color = RED;
    Material material2 = LoadMaterialDefault();
    Texture2D texture = LoadTexture("resources/textures/mass-spring.png");
    SetMaterialTexture(&material2, MATERIAL_MAP_DIFFUSE, texture);

    // material2.maps[MATERIAL_MAP_DIFFUSE].color = GREEN;
    Material material3 = LoadMaterialDefault();
    material3.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;

    RenderMesh render_mesh = RenderMesh("resources/obj/cone.obj");
    SimulationMesh sim_mesh = SimulationMesh(render_mesh);

    Mesh raylib_mesh = SimulationMesh_to_RaylibMesh(sim_mesh);

    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera, CAMERA_FREE);
        /// Keyboard controls
        if (IsKeyPressed(KEY_Q)) break;

        if (IsKeyPressed(KEY_H)) simulation_paused = !simulation_paused;
        if (IsKeyPressed(KEY_R)) state = simulation.initial_state;
        static double simulation_time = 0;
        if (!simulation_paused or IsKeyPressed(KEY_G)) {
            Clock clock = Clock(simulation_time);
            EnergyAndDerivatives f(0);
            simulation_step(simulation, state, f);
            // std::cout << "Energy " << f.energy << std::endl;
            std::cout << "Simulation time " << simulation_time << " (ms)" << std::endl;

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
                Matrix rb1_m = matrix_eigen_to_raylib(rb1.get_transformation_matrix(state));
                DrawMesh(RB_mesh, material1, rb1_m);
                simulation_render_particles(simulation.simulables, state);
                simulation_render_energies(simulation.energies, state);
                // DrawMesh(raylib_mesh, material2, MatrixIdentity());
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
