#include <algorithm>
#include <imgui.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <raylib.h>
#include <rlgl.h>
#include <vector>

#include "raylib_imgui.hpp"
#include "editor/edit_mode.hpp"
#include "linear_algebra.hpp"
#include "physics_render.hpp"
#include "raymath.h"
#include "simulable_generator.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"
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
            {
                // TODO
            }
                break;
            case FEM:
            {
                // TODO
            }
                break;
        }
    }

    enum PHYSICS_TYPE {NO_PHYISICS, MASS_SPRING, RIGID_BODY, FEM};
    PHYSICS_TYPE physics_type = NO_PHYISICS;

    private:
        PhysicsGUIGenerator phyiscs_generator;
};

struct MeshManager {
    public:
        ~MeshManager() {
            UnloadTexture(mass_spring_texture);
            UnloadTexture(rigid_body_texture);
            UnloadTexture(fem_texture);
        }
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

        void render_meshes() {
            for (const auto& pair: mesh_map) {
                switch (pair.second.physics_type) {
                    case PhysicsMesh::NO_PHYISICS:
                        DrawModel(LoadModelFromMesh(pair.second.mesh), Vector3{0}, 1, BLUE);
                        break;
                    case PhysicsMesh::MASS_SPRING:
                    {
                        Material m = LoadMaterialDefault();
                        SetMaterialTexture(&m, MATERIAL_MAP_ALBEDO, mass_spring_texture);
                        DrawMesh(pair.second.mesh, m, MatrixIdentity());
                    }
                        break;
                    case PhysicsMesh::RIGID_BODY:
                        break;
                    case PhysicsMesh::FEM:
                        break;
                }
            }
        }

    private:
        std::unordered_map<std::string, PhysicsMesh> mesh_map;
        Texture2D mass_spring_texture = LoadTexture("img/textures/mass-spring.png");
        Texture2D rigid_body_texture = LoadTexture("img/textures/rigid-body.png");
        Texture2D fem_texture = LoadTexture("img/textures/fem.png");
};

void edit_mode_sidebar(GUI_STATE& gui_state, MeshManager& mesh_manager, Simulation& out_sim, PhysicsRenderers& out_renderers) {
    gui_state = EDIT_MODE;
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport(), ImGuiDockNodeFlags_PassthruCentralNode);                // edit_mode_sidebar(gui_state, mesh_manager, simulation, physics_renderers);

    // Create side window
    ImGui::Begin("Edit scene", NULL, ImGuiWindowFlags_None);
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
        {
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
        }
            break;
        case 1: // Sphere
        {
            static int n_spheres = 0;
            if (ImGui::Button("Add")) {
                Mesh sphere = LoadMeshTinyOBJ("img/obj/sphere.obj");
                mesh_manager.add_mesh("sphere"+std::to_string(n_spheres), sphere);
                n_spheres++;
            }
        }
            break;
        case 2: // Cube
        {

            static int n_cubes = 0;
            static float sides[] = {1,1,1};
            ImGui::InputFloat3("Sizes", sides);
            if (ImGui::Button("Add")) {
                Mesh cube = GenMeshCube(sides[0], sides[1], sides[2]);
                mesh_manager.add_mesh("cube"+std::to_string(n_cubes), cube);
                n_cubes++;
            }
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
            {
                static float mass = 10;
                static float k_tension = 100;
                static float k_bending = 1;
                if (ImGui::InputFloat("mass", &mass)) mass = std::clamp(mass, 0.0f, 1000.0f);
                if (ImGui::InputFloat("tension stiffness", &k_tension)) k_tension = std::clamp(k_tension, 0.0f, 1000.0f);
                if (ImGui::InputFloat("bending stiffness", &k_bending)) k_bending = std::clamp(k_bending, 0.0f, 1000.0f);
                if (ImGui::Button("Add mass spring behaviour to mesh")) {
                    MassSpringGUIGenerator generator = {mass, k_tension, k_bending};
                    mesh_manager.add_physics_to_mesh(mesh_list[selected_mesh], generator);
                }
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
            for (const auto& pair : mesh_manager.get_mesh_map()) {
                pair.second.generate_phyiscs(out_sim, out_renderers);
                gui_state = SIMULATION_MODE;
            }
        }
    }

    ImGui::End();
}


void edit_mode_loop() {
    Camera3D camera = create_camera();
    MeshManager mesh_manager;
    Simulation simulation;
    PhysicsRenderers physics_renderers;
    GUI_STATE gui_state = EDIT_MODE;
    // Main game loop
    while (!WindowShouldClose() && gui_state == EDIT_MODE)
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera, CAMERA_ORBITAL);
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        rlDisableBackfaceCulling();
        BeginDrawing();
        {
            ClearBackground(RAYWHITE);

            BeginMode3D(camera);
            {
                DrawGrid(30, 1.0f);
                mesh_manager.render_meshes();
            }
            EndMode3D();

            DrawFPS(GetScreenWidth()*0.95, 10);
            ImGuiBeginDrawing();
            {
                edit_mode_sidebar(gui_state, mesh_manager, simulation, physics_renderers);
            }
            ImGuiEndDrawing();
        }
        EndDrawing();
        //----------------------------------------------------------------------------------
    }
    if (WindowShouldClose()) return;

    if (gui_state == SIMULATION_MODE) {
        simulation_visualization_loop(simulation, physics_renderers);
    }
}
