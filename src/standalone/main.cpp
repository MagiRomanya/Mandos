#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "render/simulation_visualization.hpp"
#include "rigid_body.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"

int main(int argc, char *argv[]) {
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos simulator");
    // Define the camera to look into our 3d world
    Camera3D camera = create_camera();

    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------

    Mesh cloth_mesh = GenMeshPlane(8.0, 4.0, 20, 20);
    const unsigned int nDoF = cloth_mesh.vertexCount * 3;
    const unsigned int n_indices = cloth_mesh.triangleCount * 3;
    std::vector<Scalar> vertices(cloth_mesh.vertices, cloth_mesh.vertices + nDoF);
    std::vector<unsigned int> indices(cloth_mesh.indices, cloth_mesh.indices + n_indices);
    // std::cout << "Number of degrees of freedom: " << nDoF << std::endl;

    Simulation simulation;
    const Scalar RB_MASS = 10.0;
    const Mat3 inertia_tensor = compute_initial_inertia_tensor(RB_MASS, indices, vertices, PARTICLES);
    simulation.initial_state.add_size(6);
    RigidBody rb(0, RB_MASS, inertia_tensor);
    simulation.simulables.rigid_bodies.push_back(rb);
    // Initial conditions
    simulation.initial_state.x.setZero();
    simulation.initial_state.v.setZero();
    simulation.initial_state.v(5) = 0.5; // add y direction angular velocity
    simulation.initial_state.v(4) = 0.01; // add y direction angular velocity

    PhysicsState state = simulation.initial_state;

    Material material = LoadMaterialDefault();
    Texture2D texture = LoadTexture("img/textures/rigid-body.png");
    SetMaterialTexture(&material, MATERIAL_MAP_ALBEDO, texture);

    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        // UpdateCamera(&camera, CAMERA_FREE);
        /// Keyboard controls
        if (IsKeyPressed(KEY_Q)) break;

        simulation_step(simulation, state);
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        rlDisableBackfaceCulling();
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);
            {
                DrawGrid(30, 1.0f);
                Mat3 rotation = rb.compute_rotation_matrix(state);
                Matrix rb_transform = { rotation(0,0), rotation(0,1), rotation(0,2), 0.0f,
                                        rotation(1,0), rotation(1,1), rotation(1,2), 0.0f,
                                        rotation(2,0), rotation(2,1), rotation(2,2), 0.0f,
                                        0.0f,          0.0f,          0.0f,          1.0f };
                DrawMesh(cloth_mesh, material, rb_transform);
            }
            EndMode3D();

            DrawFPS(10, 10);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // SimulableBounds cloth_bounds = generate_mass_spring(simulation, vertices, indices, 0.1, 100.0, 1.0);
    // MassSpringRenderer cloth_renderer = MassSpringRenderer(cloth_mesh, cloth_bounds);
    // PhysicsRenderers renderers;
    // renderers.mass_spring.push_back(cloth_renderer);

    // // froze degrees of freedom
    // simulation.frozen_dof.push_back(0);
    // simulation.frozen_dof.push_back(1);
    // simulation.frozen_dof.push_back(2);

    {
        // Test volume equation
        Mesh mesh = LoadMeshTinyOBJ("img/obj/sphere.obj");
        const unsigned int n_vertices = mesh.vertexCount * 3;
        const unsigned int n_indices = mesh.triangleCount * 3;
        const std::vector<Scalar> vertices(mesh.vertices, mesh.vertices + n_vertices);
        const std::vector<unsigned int> indices(mesh.indices, mesh.indices + n_indices);
        const Scalar volume = compute_mesh_volume(indices, vertices);
        std::cout << "Sphere volume " << volume << std::endl;
    }

    // simulation_visualization_loop(simulation, renderers);

    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
