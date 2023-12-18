#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>

#include "linear_algebra.hpp"
#include "mandos.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "physics_state.hpp"
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

    const ParticleHandle p1 = ParticleHandle(simulation, MASS)
        .add_gravity(simulation, GRAVITY)
        .set_initial_position(simulation, Vec3(0,1,0));

    const ParticleHandle p2 = ParticleHandle(simulation, MASS)
        .freeze(simulation)
        .set_initial_position(simulation, Vec3(0,1,1));

    const ParticleHandle p3 = ParticleHandle(simulation, MASS)
        .freeze(simulation)
        .set_initial_position(simulation, Vec3(0,1,-1));

    const RigidBodyHandle rb1 = RigidBodyHandle(simulation, MASS, inertia_tensor)
        .add_gravity(simulation, GRAVITY)
        .set_COM_initial_position(simulation, Vec3(0,0,0));

    ParticleRigidBodyCopuling copuling = ParticleRigidBodyCopuling(rb1.rb, p1.particle, p1.particle.get_position(simulation.initial_state.x));
    simulation.copulings.push_back(copuling);
    ParticleRigidBodyCopuling copuling2 = ParticleRigidBodyCopuling(rb1.rb, p3.particle, p3.particle.get_position(simulation.initial_state.x));
    simulation.copulings.push_back(copuling2);

    std::cout << rb1.rb.get_COM_position(simulation.initial_state.x).transpose() << std::endl;
    std::cout << rb1.rb.get_COM_position(simulation.initial_state.x_old).transpose() << std::endl;
    std::cout << simulation.initial_state.x.transpose() << std::endl;

    join_particles_with_spring(simulation, p1, p2, 10.0, 1);

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
                Matrix rb1_m = matrix_eigen_to_raylib(rb1.get_RB_transformation(state));
                DrawMesh(RB_mesh, material1, rb1_m);
                simulation_render_particles(simulation.simulables, state);
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
