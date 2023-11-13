#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <raylib.h>
#include <rlgl.h>
#include <rlImGui.h>
#include <imgui.h>

#include "linear_algebra.hpp"
#include "physics_render.hpp"
#include "physics_state.hpp"
#include "raymath.h"
#include "simulation.hpp"
#include "simulable_generator.hpp"
#include "integrators.hpp"

Camera3D create_camera(unsigned int FPS = 60) {
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

enum GUI_STATE {
EDIT_MODE,
SIMULATE_MODE,
};

int main(int argc, char *argv[]) {
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "raylib [core] example - 3d camera mode");
    rlImGuiSetup(false);
    Camera3D camera = create_camera();

    GUI_STATE gui_state = EDIT_MODE;
    std::unordered_map<std::string, Mesh> meshes;

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera, CAMERA_ORBITAL);
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        rlDisableBackfaceCulling();
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);

                DrawGrid(30, 1.0f);
                for (const auto& pair: meshes) {
                    DrawModel(LoadModelFromMesh(pair.second), Vector3{0}, 1, BLUE);
                }

            EndMode3D();

            DrawFPS(10, 10);

            rlImGuiBegin();
            {
                ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
                ImGui::SetNextWindowSizeConstraints(ImVec2(200, ImGui::GetIO().DisplaySize.y),
                                                    ImVec2(500, ImGui::GetIO().DisplaySize.y));

                // Create side window
                ImGui::Begin("Edit scene", NULL, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);
                ImGui::SeparatorText("Meshes");

                static int selected_mesh = -1;
                const char* mesh_names[] = { "Plane", "Sphere", "Cube"};

                if (ImGui::Button("Select a mesh"))
                    ImGui::OpenPopup("select_mesh_popup");
                ImGui::SameLine();
                ImGui::TextUnformatted(selected_mesh == -1 ? "<None>" : mesh_names[selected_mesh]);
                if (ImGui::BeginPopup("select_mesh_popup"))
                {
                    ImGui::SeparatorText("Available meshes");
                    for (int i = 0; i < IM_ARRAYSIZE(mesh_names); i++)
                        if (ImGui::Selectable(mesh_names[i]))
                            selected_mesh = i;
                    ImGui::EndPopup();
                }
                switch (selected_mesh) {
                    case 0: // Plane
                        static int n_divisions = 20;
                        static int n_planes = 0;
                        if (ImGui::InputInt("Plane subdivisions", &n_divisions)) n_divisions = std::clamp(n_divisions, 1, 30);
                        if (ImGui::Button("Add")) {
                            Mesh plane = GenMeshPlane(10.0f, 10.0f, n_divisions, n_divisions);
                            meshes["plane"+std::to_string(n_planes)] = plane;
                            n_planes++;
                        }
                        break;
                    case 1: // Sphere
                        static int n_spheres = 0;
                        static float radius = 1;
                        static int rings_slices[] = {20, 20};
                        if (ImGui::InputFloat("radius", &radius)) radius = std::clamp(radius, 0.1f, 10.0f);
                        ImGui::InputInt2("rings and slices", rings_slices);
                        if (ImGui::Button("Add")) {
                            Mesh sphere = GenMeshSphere(radius, 20, 20);
                            meshes["sphere"+std::to_string(n_spheres)] = sphere;
                            n_spheres++;
                        }
                        break;
                    case 2: // Cube
                        static int n_cubes = 0;
                        static float sides[] = {1,1,1};
                        ImGui::InputFloat3("Sizes", sides);
                        if (ImGui::Button("Add")) {
                            Mesh cube = GenMeshCube(sides[0], sides[1], sides[2]);
                            meshes["cube"+std::to_string(n_spheres)] = cube;
                            n_spheres++;
                        }
                        break;
                }
                const char* mesh_list[20];
                int i = 0;
                for (const auto& pair : meshes) {
                    mesh_list[i++] = pair.first.c_str();
                }
                static int mesh_current = -1;
                ImGui::ListBox("Meshes", &mesh_current, mesh_list, meshes.size(), 4);

                if (mesh_current >= 0) {
                    if (ImGui::Button("Delete")) {
                        meshes.erase(mesh_list[mesh_current]);
                        mesh_current = -1;
                    }
                    static int e = 0;
                    ImGui::RadioButton("Mass Spring", &e, 0); ImGui::SameLine();
                    ImGui::RadioButton("Rigid Body", &e, 1); ImGui::SameLine();
                    ImGui::RadioButton("FEM", &e, 2);
                    switch (e) {
                        case 0: // Mass Spring
                            static float mass = 1;
                            static float k_tension = 100;
                            static float k_bending = 1;
                            if (ImGui::InputFloat("mass", &mass)) mass = std::clamp(mass, 0.0f, 1000.0f);
                            if (ImGui::InputFloat("tension stiffness", &k_tension)) k_tension = std::clamp(k_tension, 0.0f, 1000.0f);
                            if (ImGui::InputFloat("bending stiffness", &k_bending)) k_bending = std::clamp(k_bending, 0.0f, 1000.0f);
                            break;
                        case 1: // Rigid Body
                            ImGui::Text("TODO!");
                            break;
                        case 2: // FEM
                            ImGui::Text("TODO!");
                            break;
                    }
                }

                if (ImGui::TreeNode("Tree Node"))
                {
                    ImGui::Text("Tree node content");
                    ImGui::TreePop();
                }
                ImGui::End();
                bool open = true;
                ImGui::ShowDemoWindow(&open);

            }
            rlImGuiEnd();
        EndDrawing();
        //----------------------------------------------------------------------------------
    }
    // De-Initialization
    //--------------------------------------------------------------------------------------
    rlImGuiShutdown();
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
