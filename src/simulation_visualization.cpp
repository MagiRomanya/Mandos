#include "render/simulation_visualization.hpp"
#include "imgui.h"
#include "physics_state.hpp"
#include "raylib.h"
#include <iostream>

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

    bool simulation_pause = false;
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

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);

                DrawGrid(30, 1.0f);
                phy_renderers.draw(state);

            EndMode3D();

            DrawFPS(10, 10);

            rlImGuiBegin();
            {
                // show ImGui Content
                const float screen_width = ImGui::GetIO().DisplaySize.x;
                const float screen_hieght = ImGui::GetIO().DisplaySize.y;
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 200, 0), ImGuiCond_Always);

    // Set the window size
    ImGui::SetNextWindowSize(ImVec2(200, screen_hieght), ImGuiCond_Always);

    // Begin the window
    ImGui::Begin("Right-aligned Window", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
                // ImGui::SetNextWindowSizeConstraints(ImVec2(200, screen_hieght),
                //                                     ImVec2(500, screen_hieght));

                // Create side window
                // ImGui::Begin("Simulation manager", NULL, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);

                ImGui::End();
                ImGui::ShowDemoWindow();
            }
            rlImGuiEnd();
        EndDrawing();
        //----------------------------------------------------------------------------------
    }
}
