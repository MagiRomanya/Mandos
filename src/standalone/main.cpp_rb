#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>
#include <Eigen/LU>


#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "raymath.h"
#include "render/simulation_visualization.hpp"
#include "rigid_body.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"
#include "render/draw_vector.hpp"

int main(int argc, char *argv[]) {
    FEM_ElementParameters delete_this_variable(1,1,Vec4::Zero(), Eigen::Matrix<Scalar,4,3>::Zero());
    FEM_Element3D elemenr(Particle(1,1), Particle(1,1), Particle(1,1), Particle(1,1), delete_this_variable);

    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos simulator");
    // Define the camera to look into our 3d world
    Camera3D camera = { 0 };
    camera.position = Vector3{ 0.0f, 10.f, 10.0f };  // Camera position
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type

    // Make the window resizable
    SetWindowState(FLAG_WINDOW_RESIZABLE);

    SetTargetFPS(200);               // Set our game to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------

    Mesh cloth_mesh = GenMeshPlane(8.0, 1.0, 20, 20);
    Mesh vector_mesh = LoadMeshTinyOBJ("img/obj/vector.obj");
    // Mesh cloth_mesh = LoadMeshTinyOBJ("img/obj/hammer.obj");
    const unsigned int nDoF = cloth_mesh.vertexCount * 3;
    const unsigned int n_indices = cloth_mesh.triangleCount * 3;
    std::vector<Scalar> vertices(cloth_mesh.vertices, cloth_mesh.vertices + nDoF);
    std::vector<unsigned int> indices(cloth_mesh.indices, cloth_mesh.indices + n_indices);
    // std::cout << "Number of degrees of freedom: " << nDoF << std::endl;

    Simulation simulation;
    const Scalar RB_MASS = 1.0;
    const Mat3 inertia_tensor = compute_initial_inertia_tensor(RB_MASS, indices, vertices, PARTICLES);
    simulation.initial_state.add_size(6);
    RigidBody rb(0, RB_MASS, inertia_tensor);
    simulation.simulables.rigid_bodies.push_back(rb);
    // Initial conditions
    simulation.initial_state.x.setZero();
    simulation.initial_state.v.setZero();
    simulation.initial_state.v(5) = 1; // add y direction angular velocity
    simulation.initial_state.v(4) = 0.001; // add y direction angular velocity

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
                // std::cout << rotation.determinant() << std::endl;
                // std::cout << rotation.transpose() * rotation << std::endl;
                Matrix rb_transform = { rotation(0,0), rotation(0,1), rotation(0,2), 0.0f,
                                        rotation(1,0), rotation(1,1), rotation(1,2), 0.0f,
                                        rotation(2,0), rotation(2,1), rotation(2,2), 0.0f,
                                        0.0f,          0.0f,          0.0f,          1.0f };

                Matrix vector_transform = MatrixMultiply(MatrixRotateX(3.14159 / 2), rb_transform);
                DrawMesh(cloth_mesh, material, rb_transform);
                DrawMesh(vector_mesh, material, vector_transform);

                const Vec3 theta = state.x.tail(3);
                const Vec3 omega = state.v.tail(3);
                DrawVector(Vec3::Zero(), theta, BLUE, theta.norm()/ 10);
                DrawVector(Vec3::Zero(), omega, RED);
                // DrawVector(Vec3::Zero(), Vec3(1,0,0), RED);
                // DrawVector(Vec3::Zero(), Vec3(0,1,0), GREEN);
                // DrawVector(Vec3::Zero(), Vec3(0,0,1), BLUE);
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
