#include <imgui.h>
#include <raylib.h>
#include <rlgl.h>
#include <iostream>

#include "raylib_imgui.hpp"
#include "render/simulation_visualization.hpp"
#include "physics_state.hpp"

Camera3D create_camera(unsigned int FPS) {
    // Define the camera to look into our 3d world
    Camera3D camera = { 0 };
    camera.position = Vector3{ 0.0f, 4.0f, 10.0f };  // Camera position
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type

    SetTargetFPS(FPS);               // Set our game to run at 60 frames-per-second

    // Make the window resizable
    SetWindowState(FLAG_WINDOW_RESIZABLE);
    return camera;
}

void simulation_visualization_loop(Simulation& simulation, PhysicsRenderers& phy_renderers) {
    Camera3D camera = create_camera();

    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second

    PhysicsState state = simulation.initial_state;

    bool simulation_pause = true;
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera, CAMERA_ORBITAL);
        /// Keyboard controls
        if (IsKeyPressed(KEY_Q)) break;

        EnergyAndDerivatives f(0);
        if (!simulation_pause)
            simulation_step(simulation, state, f);
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        rlDisableBackfaceCulling();
        BeginDrawing();
        {
            ClearBackground(RAYWHITE);

            BeginMode3D(camera);

                DrawGrid(30, 1.0f);
                phy_renderers.draw(state);

            EndMode3D();

            DrawFPS(10, 10);

            ImGuiBeginDrawing();
            {
                // show ImGui Content
                const float screen_width = ImGui::GetIO().DisplaySize.x;
                const float screen_hieght = ImGui::GetIO().DisplaySize.y;

                // Create side window
                ImGui::DockSpaceOverViewport(ImGui::GetMainViewport(), ImGuiDockNodeFlags_PassthruCentralNode);                // edit_mode_sidebar(gui_state, mesh_manager, simulation, physics_renderers);
                ImGui::Begin("Simulation manager", NULL, ImGuiWindowFlags_None);

                ImGui::End();
                ImGui::ShowDemoWindow();
            }
            ImGuiEndDrawing();
;        }
        EndDrawing();
        //----------------------------------------------------------------------------------
    }
}
