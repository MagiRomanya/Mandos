#ifndef SIMULATION_VISUALIZATION_H_
#define SIMULATION_VISUALIZATION_H_

#include <rlImGui.h>
#include <rlgl.h>
#include <imgui.h>

#include "physics_render.hpp"
#include "physics_state.hpp"
#include "simulation.hpp"

void simulation_visualization_loop(Simulation& simulation, PhysicsRenderers& phy_renderers) {
    Camera3D camera = { 0 };
    camera.position = Vector3{ 0.0f, 4.0f, 10.0f };  // Camera position
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type

    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second

    PhysicsState state = simulation.initial_state;

    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera, CAMERA_ORBITAL);
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

                DrawGrid(10, 1.0f);
                phy_renderers.draw(state);

            EndMode3D();

            DrawFPS(10, 10);

            rlImGuiBegin();
            {
                // show ImGui Content
                // bool open = true;
                // ImGui::ShowDemoWindow(&open);
            }
            rlImGuiEnd();
        EndDrawing();
        //----------------------------------------------------------------------------------
    }
}

#endif // SIMULATION_VISUALIZATION_H_
