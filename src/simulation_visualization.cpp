#include <imgui.h>
#include <raylib.h>
#include <rlgl.h>
#include <iostream>
#include <vector>

#include "linear_algebra.hpp"
#include "raylib_imgui.hpp"
#include "render/simulation_visualization.hpp"
#include "physics_state.hpp"
#include "clock.hpp"
#include "rigid_body.hpp"
#include "spring.hpp"

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

void simulation_render_simulables_and_energies(const Simulation& simulation, const PhysicsState& state) {
    simulation_render_energies(simulation.energies, state);
    simulation_render_simulables(simulation.simulables, state);
}

void simulation_render_simulables(const Simulables& simulables, const PhysicsState& state) {
    // DRAW PARTICLES
    //----------------------------------------------------------------------------------
    Color colors[] = {LIGHTGRAY, GRAY, DARKGRAY, YELLOW, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, SKYBLUE, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BEIGE, BROWN, DARKBROWN};
    for (unsigned int i = 0; i < simulables.particles.size(); i++) {
        const Particle& p = simulables.particles[i];
        const Vec3 x = p.get_position(state);
        DrawSphere(Vector3{x.x(), x.y(), x.z()}, 0.05, colors[i % IM_ARRAYSIZE(colors)]);
    }

    // DRAW RIGID BODIES
    //----------------------------------------------------------------------------------
    for (unsigned int i = 0; i < simulables.rigid_bodies.size(); i++) {
        const RigidBody& rb = simulables.rigid_bodies[i];
        const Vec3 x = rb.get_COM_position(state);
        const Mat3 rot = rb.compute_rotation_matrix(state);
        std::cout << "DRAW RIGID BODIES NOT IMPLEMENTED YET" << std::endl;
    }
}
void simulation_render_energies(const Energies& energies, const PhysicsState& state) {
    // DRAW SPRINGS
    //----------------------------------------------------------------------------------
    for (unsigned int i = 0; i < energies.particle_springs.size(); i++) {
        const ParticleSpring& s = energies.particle_springs[i];
        const Vec3& x1 = s.p1.get_position(state);
        const Vec3& x2 = s.p2.get_position(state);
        DrawLine3D(Vector3{x1.x(),x1.y(),x1.z()}, Vector3{x2.x(),x2.y(),x2.z()}, RED);
    }
    // DRAW FEM TETRAHEDRONS
    //----------------------------------------------------------------------------------
    for (unsigned int i = 0; i < energies.fem_elements_3d.size(); i++) {
        const FEM_Element3D& e = energies.fem_elements_3d[i];
        const Vec3& x1 = e.p1.get_position(state);
        const Vec3& x2 = e.p2.get_position(state);
        const Vec3& x3 = e.p3.get_position(state);
        const Vec3& x4 = e.p4.get_position(state);
        DrawLine3D(Vector3{x1.x(), x1.y(), x1.z()}, Vector3{x2.x(), x2.y(), x2.z()}, BLUE);
        DrawLine3D(Vector3{x1.x(), x1.y(), x1.z()}, Vector3{x3.x(), x3.y(), x3.z()}, BLUE);
        DrawLine3D(Vector3{x1.x(), x1.y(), x1.z()}, Vector3{x4.x(), x4.y(), x4.z()}, BLUE);
        DrawLine3D(Vector3{x4.x(), x4.y(), x4.z()}, Vector3{x2.x(), x2.y(), x2.z()}, BLUE);
        DrawLine3D(Vector3{x4.x(), x4.y(), x4.z()}, Vector3{x3.x(), x3.y(), x3.z()}, BLUE);
        DrawLine3D(Vector3{x2.x(), x2.y(), x2.z()}, Vector3{x3.x(), x3.y(), x3.z()}, BLUE);
        const Vec3 com = (x1 +x2 +x3 +x4) / 4;
        DrawSphere(Vector3{com.x(), com.y(),com.z()}, 0.04, PINK);
    }
}
