#include <imgui.h>
#include <raylib.h>
#include <rlgl.h>
#include <iostream>
#include <vector>

#include "raylib_imgui.hpp"
#include "render/simulation_visualization.hpp"
#include "physics_state.hpp"
#include "clock.hpp"

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

enum SIMULATION_FLOW_STATE { SIMULATION_PAUSE, SIMULATION_STEP, SIMULATION_RUN };

void simulation_visualization_loop(Simulation& simulation, PhysicsRenderers& phy_renderers) {
    Camera3D camera = create_camera();

    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second

    PhysicsState state = simulation.initial_state;

    SIMULATION_FLOW_STATE simulation_flow_state = SIMULATION_PAUSE;
    unsigned int n_simulation_steps = 0;

    EnergyAndDerivatives f(0);
    bool exit_simulation_enviroment = false;
    double simulation_time = 0;
    while (!WindowShouldClose() && !exit_simulation_enviroment)    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        // UpdateCamera(&camera, CAMERA_ORBITAL);

        /// Keyboard controls
        if (IsKeyPressed(KEY_Q)) break;

        // SIMULATION FLOW
        switch (simulation_flow_state) {
        case SIMULATION_PAUSE:
            break;
        case SIMULATION_RUN:
            {
                Clock clock(simulation_time);
                simulation_step(simulation, state, f);
                n_simulation_steps++;
            }
            break;
        case SIMULATION_STEP:
            {
                Clock clock(simulation_time);
                simulation_step(simulation, state, f);
                n_simulation_steps++;
                simulation_flow_state = SIMULATION_PAUSE;
            }
            break;
        }
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
                // Create side window
                ImGui::DockSpaceOverViewport(ImGui::GetMainViewport(), ImGuiDockNodeFlags_PassthruCentralNode);
                ImGui::Begin("Simulation manager", NULL, ImGuiWindowFlags_None);

                ImGui::SeparatorText("Simulation control");
                ImGui::Text("Simulation steps = %u", n_simulation_steps);
                {
                    static bool simulation_run = false;
                    if (ImGui::Button("Step")) {
                        simulation_flow_state = SIMULATION_STEP;
                        simulation_run = false;
                    }
                    if (ImGui::Checkbox("Simulation run", &simulation_run)) {
                        if (simulation_run) simulation_flow_state = SIMULATION_RUN;
                        else simulation_flow_state = SIMULATION_PAUSE;
                    }
                }

                if (ImGui::Button("Reset to initial state")) {
                    state = simulation.initial_state;
                    n_simulation_steps = 0;
                }

                ImGui::SeparatorText("Analytics");
                static std::vector<float> energy_list;
                static std::vector<float> simulation_time_list;

                if (ImGui::Button("Discard all data")) {
                    energy_list.clear();
                    simulation_time_list.clear();
                }
                if ((simulation_flow_state == SIMULATION_RUN) or
                    (simulation_flow_state == SIMULATION_STEP)) {
                    const float energy = f.energy;
                    energy_list.push_back(energy);
                    simulation_time_list.push_back(simulation_time);
                }
                ImGui::PlotLines("Total energy", energy_list.data(), energy_list.size(), 0, NULL, FLT_MAX, FLT_MAX, ImVec2(0.0, 80.0));
                if (ImGui::Button("Discard energy data")) energy_list.clear();

                ImGui::PlotLines("Step time cost (ms)", simulation_time_list.data(), simulation_time_list.size(), 0, NULL, FLT_MAX, FLT_MAX, ImVec2(0.0, 80.0));
                if (ImGui::Button("Discard time cost data")) simulation_time_list.clear();

                if (ImGui::Button("Exit")) {
                    exit_simulation_enviroment = true;
                }
                ImGui::End();

            }
            ImGuiEndDrawing();
;        }
        EndDrawing();
        //----------------------------------------------------------------------------------
    }
}
