#include <cstring>
#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "render/simulation_visualization.hpp"
#include "simulable_generator.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"

int main(void) {
    // Tetgen test
    // std::vector<float> vertices;
    // std::vector<unsigned int> indices;
    // LoadVerticesAndIndicesTinyOBJ("img/obj/cube.obj", vertices, indices);


    // std::vector<float> out_vertex;
    // std::vector<unsigned int> out_index;
    // tetgen_compute_tetrahedrons(indices, vertices, out_index, out_vertex);
    // Simulation simulation;
    // const Scalar young_modulus = 50e3;
    // const Scalar poisson_ratio = 0.2;
    // // generate_FEM3D_tetrahedron(simulation, 1, poisson_ratio, young_modulus);
    // generate_FEM3D_from_tetrahedron_mesh(simulation, 1, poisson_ratio, young_modulus, out_index, out_vertex);

    // // Initialization
    // //--------------------------------------------------------------------------------------
    // const int screenWidth = 1600;
    // const int screenHeight = 900;

    // InitWindow(screenWidth, screenHeight, "Mandos simulator");
    // Camera3D camera = create_camera();

    // SetTargetFPS(200);
    // //--------------------------------------------------------------------------------------

    // bool simulation_paused = true;
    // PhysicsState state = simulation.initial_state;
    // while (!WindowShouldClose())    // Detect window close button or ESC key
    // {
    //     // Update
    //     //----------------------------------------------------------------------------------
    //     UpdateCamera(&camera, CAMERA_ORBITAL);
    //     /// Keyboard controls
    //     if (IsKeyPressed(KEY_Q)) break;

    //     if (IsKeyPressed(KEY_H)) simulation_paused = !simulation_paused;
    //     if (IsKeyPressed(KEY_R)) state = simulation.initial_state;
    //     if (!simulation_paused or IsKeyPressed(KEY_G)) {
    //         EnergyAndDerivatives f(0);
    //         simulation_step(simulation, state, f);
    //         std::cout << "Energy " << f.energy << std::endl;
    //     }
    //     //----------------------------------------------------------------------------------

    //     // Draw
    //     //----------------------------------------------------------------------------------
    //     BeginDrawing();

    //         ClearBackground(RAYWHITE);

    //         BeginMode3D(camera);
    //         {
    //             DrawGrid(30, 1.0f);
    //             // DrawPoint3D(Vector3{0, 1, 0}, RED);
    //             // for (unsigned int i = 0; i < out_vertex.size()/3; i++) {
    //             //     DrawPoint3D(Vector3{out_vertex[3*i], out_vertex[3*i+1], out_vertex[3*i+2]}, RED);
    //             // }
    //             // for (unsigned int i = 0; i < out_index.size()/4; i++) {
    //             //     const Vector3 p1 = Vector3{out_vertex[3*out_index[4*i+0]], out_vertex[3*out_index[4*i+0]+1], out_vertex[3*out_index[4*i+0]+2]};
    //             //     const Vector3 p2 = Vector3{out_vertex[3*out_index[4*i+1]], out_vertex[3*out_index[4*i+1]+1], out_vertex[3*out_index[4*i+1]+2]};
    //             //     const Vector3 p3 = Vector3{out_vertex[3*out_index[4*i+2]], out_vertex[3*out_index[4*i+2]+1], out_vertex[3*out_index[4*i+2]+2]};
    //             //     const Vector3 p4 = Vector3{out_vertex[3*out_index[4*i+3]], out_vertex[3*out_index[4*i+3]+1], out_vertex[3*out_index[4*i+3]+2]};
    //             //     DrawLine3D(p1, p2, BLUE);
    //             //     DrawLine3D(p1, p3, BLUE);
    //             //     DrawLine3D(p1, p4, BLUE);
    //             //     DrawLine3D(p2, p3, BLUE);
    //             //     DrawLine3D(p2, p4, BLUE);
    //             //     DrawLine3D(p3, p4, BLUE);
    //             // }
    //             simulation_render_simulables_and_energies(simulation, state);
    //         }
    //         EndMode3D();

    //         DrawFPS(10, 10);

    //     EndDrawing();
    //     //----------------------------------------------------------------------------------
    // }
    // // De-Initialization
    // //--------------------------------------------------------------------------------------
    // CloseWindow();        // Close window and OpenGL context
    // //--------------------------------------------------------------------------------------

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
                // DrawVector(Vec3::Zero(), theta, BLUE, theta.norm()/ 10);
                // DrawVector(Vec3::Zero(), omega, RED);
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
    return 0;
}
