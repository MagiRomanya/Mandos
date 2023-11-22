#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>


#include "fem_unit.hpp"
#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "render/simulation_visualization.hpp"
#include "simulable_generator.hpp"
#include "simulation.hpp"

int main(void) {
    Simulation simulation;
    const Scalar young_modulus = 50e3;
    const Scalar poisson_ratio = 0.2;
    generate_FEM3D_tetrahedron(simulation, 1, poisson_ratio, young_modulus);

    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos simulator");
    Camera3D camera = create_camera();

    SetTargetFPS(200);
    //--------------------------------------------------------------------------------------

    PhysicsState state = simulation.initial_state;
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera, CAMERA_ORBITAL);
        /// Keyboard controls
        if (IsKeyPressed(KEY_Q)) break;

        if (true or IsKeyPressed(KEY_G)) {
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
