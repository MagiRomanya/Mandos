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
    std::vector<float> vertices;
    std::vector<unsigned int> indices;
    LoadVerticesAndIndicesTinyOBJ("img/obj/cube.obj", vertices, indices);


    std::vector<float> out_vertex;
    std::vector<unsigned int> out_index;
    tetgen_compute_tetrahedrons(indices, vertices, out_index, out_vertex);
    Simulation simulation;
    const Scalar young_modulus = 50e3;
    const Scalar poisson_ratio = 0.2;
    // generate_FEM3D_tetrahedron(simulation, 1, poisson_ratio, young_modulus);
    generate_FEM3D_from_tetrahedron_mesh(simulation, 1, poisson_ratio, young_modulus, out_index, out_vertex);

    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos simulator");
    Camera3D camera = create_camera();

    SetTargetFPS(200);
    //--------------------------------------------------------------------------------------

    bool simulation_paused = true;
    PhysicsState state = simulation.initial_state;
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
                // DrawPoint3D(Vector3{0, 1, 0}, RED);
                // for (unsigned int i = 0; i < out_vertex.size()/3; i++) {
                //     DrawPoint3D(Vector3{out_vertex[3*i], out_vertex[3*i+1], out_vertex[3*i+2]}, RED);
                // }
                // for (unsigned int i = 0; i < out_index.size()/4; i++) {
                //     const Vector3 p1 = Vector3{out_vertex[3*out_index[4*i+0]], out_vertex[3*out_index[4*i+0]+1], out_vertex[3*out_index[4*i+0]+2]};
                //     const Vector3 p2 = Vector3{out_vertex[3*out_index[4*i+1]], out_vertex[3*out_index[4*i+1]+1], out_vertex[3*out_index[4*i+1]+2]};
                //     const Vector3 p3 = Vector3{out_vertex[3*out_index[4*i+2]], out_vertex[3*out_index[4*i+2]+1], out_vertex[3*out_index[4*i+2]+2]};
                //     const Vector3 p4 = Vector3{out_vertex[3*out_index[4*i+3]], out_vertex[3*out_index[4*i+3]+1], out_vertex[3*out_index[4*i+3]+2]};
                //     DrawLine3D(p1, p2, BLUE);
                //     DrawLine3D(p1, p3, BLUE);
                //     DrawLine3D(p1, p4, BLUE);
                //     DrawLine3D(p2, p3, BLUE);
                //     DrawLine3D(p2, p4, BLUE);
                //     DrawLine3D(p3, p4, BLUE);
                // }
                simulation_render_simulables_and_energies(simulation, state);
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
