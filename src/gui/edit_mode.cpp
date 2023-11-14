#include <algorithm>
#include <imgui.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <raylib.h>
#include <rlgl.h>
#include <vector>

#include "gui/edit_mode.hpp"
#include "linear_algebra.hpp"
#include "physics_render.hpp"
#include "rlImGui.h"
#include "simulable_generator.hpp"
#include "simulation.hpp"
#include "render/simulation_visualization.hpp"

struct PhysicsMesh {
    Mesh mesh;
    void add_physics(MassSpringGUIGenerator gen) {
        physics_type = MASS_SPRING;
        phyiscs_generator.mass_spring = gen;
    }

    void remove_physics() { physics_type = NO_PHYISICS; }

    void generate_phyiscs(Simulation& simulation, PhysicsRenderers& physics_renderers) const {
        const unsigned int nDoF = mesh.vertexCount * 3;
        const unsigned int n_indices = mesh.triangleCount * 3;
        std::vector<Scalar> vertices(mesh.vertices, mesh.vertices + nDoF);
        std::vector<unsigned int> indices(mesh.indices, mesh.indices + n_indices);

        switch (physics_type) {
            case NO_PHYISICS:
                break;
            case MASS_SPRING:
            {
                /// GENERATE PHYISICS
                const MassSpringGUIGenerator& g = phyiscs_generator.mass_spring;
                SimulableBounds bounds = generate_mass_spring(simulation, vertices, indices, g.mass *3 / nDoF, g.k_tension, g.k_bending);
                // Add the frozen nodes with the propper dof offset
                for (unsigned int i = 0; i < g.frozen_nodes.size(); i++) {
                    simulation.frozen_dof.push_back(g.frozen_nodes[i] + bounds.dof_index);
                }

                /// GENERATE RENDERING
                physics_renderers.mass_spring.emplace_back(mesh, bounds);
            }
                break;
            case RIGID_BODY:
                break;
            case FEM:
                break;
        }
    }

    private:
        PhysicsGUIGenerator phyiscs_generator;
        enum PHYSICS_TYPE {NO_PHYISICS, MASS_SPRING, RIGID_BODY, FEM};
        PHYSICS_TYPE physics_type = NO_PHYISICS;
};

struct MeshManager {
    public:
        void add_mesh(std::string name, Mesh mesh) { mesh_map[name].mesh = mesh; }
        void delete_mesh(std::string name) {
            UnloadMesh(mesh_map[name].mesh);
            mesh_map.erase(name);
        }
        const std::unordered_map<std::string, PhysicsMesh>& get_mesh_map() const { return mesh_map; }
        int size() const { return mesh_map.size(); }

        void add_physics_to_mesh(std::string name, MassSpringGUIGenerator generator) {
            mesh_map[name].add_physics(generator);
        }
    private:
        std::unordered_map<std::string, PhysicsMesh> mesh_map;
};

MeshManager mesh_manager;

void edit_mode_render_meshes() {
    for (const auto& pair: mesh_manager.get_mesh_map()) {
        DrawModel(LoadModelFromMesh(pair.second.mesh), Vector3{0}, 1, BLUE);
    }
}

void edit_mode_sidebar() {
    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSizeConstraints(ImVec2(200, ImGui::GetIO().DisplaySize.y),
                                        ImVec2(500, ImGui::GetIO().DisplaySize.y));

    // Create side window
    ImGui::Begin("Edit scene", NULL, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);
    ImGui::SeparatorText("Meshes");

    static int selected_mesh_generator = -1;
    const char* mesh_names[] = { "Plane", "Sphere", "Cube"};

    if (ImGui::Button("Select a mesh"))
        ImGui::OpenPopup("select_mesh_popup");
    ImGui::SameLine();
    ImGui::TextUnformatted(selected_mesh_generator == -1 ? "<None>" : mesh_names[selected_mesh_generator]);
    if (ImGui::BeginPopup("select_mesh_popup"))
    {
        ImGui::SeparatorText("Available meshes");
        for (int i = 0; i < IM_ARRAYSIZE(mesh_names); i++)
            if (ImGui::Selectable(mesh_names[i]))
                selected_mesh_generator = i;
        ImGui::EndPopup();
    }
    static int selected_mesh = -1;
    switch (selected_mesh_generator) {
        case 0: // Plane
            static int n_divisions = 20;
            static float plane_dimensions[] = {10.0f, 10.0f};
            static int n_planes = 0;
            ImGui::InputFloat2("Plane width and height", plane_dimensions);
            if (ImGui::InputInt("Plane subdivisions", &n_divisions)) n_divisions = std::clamp(n_divisions, 1, 30);
            if (ImGui::Button("Add")) {
                Mesh plane = GenMeshPlane(10.0f, 10.0f, n_divisions, n_divisions);
                mesh_manager.add_mesh("plane"+std::to_string(n_planes), plane);
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
                mesh_manager.add_mesh("sphere"+std::to_string(n_spheres), sphere);
                n_spheres++;
            }
            break;
        case 2: // Cube
            static int n_cubes = 0;
            static float sides[] = {1,1,1};
            ImGui::InputFloat3("Sizes", sides);
            if (ImGui::Button("Add")) {
                Mesh cube = GenMeshCube(sides[0], sides[1], sides[2]);
                mesh_manager.add_mesh("cube"+std::to_string(n_cubes), cube);
                n_cubes++;
            }
            break;
    }
    const char* mesh_list[20];
    int i = 0;
    for (const auto& pair : mesh_manager.get_mesh_map()) {
        mesh_list[i++] = pair.first.c_str();
    }
    ImGui::ListBox("Meshes", &selected_mesh, mesh_list, mesh_manager.size(), 4);

    if (selected_mesh >= 0) {
        if (ImGui::Button("Delete")) {
            mesh_manager.delete_mesh(mesh_list[selected_mesh]);
            if (mesh_manager.size() == 0) selected_mesh = -1;
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
                if (ImGui::Button("Add mass spring behaviour to mesh")) {
                    MassSpringGUIGenerator generator = {mass, k_tension, k_bending};
                    mesh_manager.add_physics_to_mesh(mesh_list[selected_mesh], generator);
                }
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
    if (mesh_manager.size() != 0) {
        if (ImGui::Button("Simulate!")) {
            Simulation simulation;
            PhysicsRenderers physics_renderers;
            for (const auto& pair : mesh_manager.get_mesh_map()) {
                pair.second.generate_phyiscs(simulation, physics_renderers);
            }
            simulation_visualization_loop(simulation, physics_renderers);
        }
    }

    ImGui::End();
}

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

void edit_mode_loop() {
    Camera3D camera = create_camera();
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
                edit_mode_render_meshes();

            EndMode3D();

            DrawFPS(GetScreenWidth()*0.95, 10);

            rlImGuiBegin();
            {
                edit_mode_sidebar();
                // bool open = true;
                // ImGui::ShowDemoWindow(&open);
            }
            rlImGuiEnd();
        EndDrawing();
        //----------------------------------------------------------------------------------
    }
}
