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
    ~PhysicsMesh() {
        UnloadMesh(mesh);
    }
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
                SimulableBounds bounds = generate_mass_spring(simulation, vertices, indices, g.mass *3 / nDoF, g.k_tension, g.k_bending, g.damping);
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
                        DrawModel(LoadModelFromMesh(pair.second.mesh), Vector3{0,0,0}, 1, BLUE);
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
        Texture2D mass_spring_texture = LoadTexture("resources/textures/mass-spring.png");
        Texture2D rigid_body_texture = LoadTexture("resources/textures/rigid-body.png");
        Texture2D fem_texture = LoadTexture("resources/textures/fem.png");
};

struct EditModeUserSelectionState {
    // Mesh
    int selected_mesh_generator = -1;
    int selected_mesh = -1;

    // Geometries
    int n_divisions = 20;
    float plane_dimensions[2] = {10.0f, 10.0f};
    int n_planes = 0;
    int n_spheres = 0;
    int n_cubes = 0;

    // Simulable selector
    int simulable = 0;

    // Mass spring
    float mass = 100;
    float k_tension = 100;
    float k_bending = 1;
    float damping = 1;

    float TimeStep = 0.1;
};

void edit_mode_sidebar(EditModeUserSelectionState& s_state, GUI_STATE& gui_state, MeshManager& mesh_manager, Simulation& out_sim, PhysicsRenderers& out_renderers) {
    gui_state = EDIT_MODE;
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport(), ImGuiDockNodeFlags_PassthruCentralNode);                // edit_mode_sidebar(gui_state, mesh_manager, simulation, physics_renderers);

    // Create side window
    ImGui::Begin("Edit scene", NULL, ImGuiWindowFlags_None);
    ImGui::SeparatorText("Meshes");

    const char* mesh_names[] = { "Plane", "Sphere", "Cube"};

    if (ImGui::Button("Select a mesh"))
        ImGui::OpenPopup("select_mesh_popup");
    ImGui::SameLine();
    ImGui::TextUnformatted(s_state.selected_mesh_generator == -1 ? "<None>" : mesh_names[s_state.selected_mesh_generator]);
    if (ImGui::BeginPopup("select_mesh_popup"))
    {
        ImGui::SeparatorText("Available meshes");
        for (int i = 0; i < IM_ARRAYSIZE(mesh_names); i++)
            if (ImGui::Selectable(mesh_names[i]))
                s_state.selected_mesh_generator = i;
        ImGui::EndPopup();
    }
    switch (s_state.selected_mesh_generator) {
        case 0: // Plane
        {
            ImGui::InputFloat2("Plane width and height", s_state.plane_dimensions);
            if (ImGui::InputInt("Plane subdivisions", &s_state.n_divisions)) s_state.n_divisions = std::clamp(s_state.n_divisions, 1, 30);
            if (ImGui::Button("Add")) {
                Mesh plane = GenMeshPlane(10.0f, 10.0f, s_state.n_divisions, s_state.n_divisions);
                mesh_manager.add_mesh("plane"+std::to_string(s_state.n_planes), plane);
                s_state.n_planes++;
            }
        }
            break;
        case 1: // Sphere
        {
            if (ImGui::Button("Add")) {
                Mesh sphere = LoadMeshTinyOBJ("resources/obj/sphere.obj");
                mesh_manager.add_mesh("sphere"+std::to_string(s_state.n_spheres), sphere);
                s_state.n_spheres++;
            }
        }
            break;
        case 2: // Cube
        {

            ImGui::Text("TODO cube");
        }
            break;
    }
    const char* mesh_list[20];
    int i = 0;
    for (const auto& pair : mesh_manager.get_mesh_map()) {
        mesh_list[i++] = pair.first.c_str();
    }
    ImGui::ListBox("Meshes", &s_state.selected_mesh, mesh_list, mesh_manager.size(), 4);

    if (s_state.selected_mesh >= 0) {
        if (ImGui::Button("Delete")) {
            mesh_manager.delete_mesh(mesh_list[s_state.selected_mesh]);
            if (mesh_manager.size() == 0) s_state.selected_mesh = -1;
        }
        ImGui::RadioButton("Mass Spring", &s_state.simulable, 0); ImGui::SameLine();
        ImGui::RadioButton("Rigid Body", &s_state.simulable, 1); ImGui::SameLine();
        ImGui::RadioButton("FEM", &s_state.simulable, 2);
        switch (s_state.simulable) {
            case 0: // Mass Spring
            {
                if (ImGui::InputFloat("mass", &s_state.mass)) s_state.mass = std::clamp(s_state.mass, 0.0f, 1000.0f);
                if (ImGui::InputFloat("tension stiffness", &s_state.k_tension)) s_state.k_tension = std::clamp(s_state.k_tension, 0.0f, 1000.0f);
                if (ImGui::InputFloat("bending stiffness", &s_state.k_bending)) s_state.k_bending = std::clamp(s_state.k_bending, 0.0f, 1000.0f);
                if (ImGui::InputFloat("Damping", &s_state.damping)) s_state.damping = std::clamp(s_state.damping, 0.0f, 1000.0f);
                if (ImGui::Button("Add mass spring behaviour to mesh")) {
                    MassSpringGUIGenerator generator = {s_state.mass, s_state.k_tension, s_state.k_bending, s_state.damping};
                    mesh_manager.add_physics_to_mesh(mesh_list[s_state.selected_mesh], generator);
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
    if (mesh_manager.size() != 0) {
        ImGui::SeparatorText("Simulation settings");
        ImGui::InputFloat("Delta Time", &s_state.TimeStep);
        if (ImGui::Button("Simulate!")) {
            out_sim.TimeStep = s_state.TimeStep;
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
    EditModeUserSelectionState selection_state;
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
                edit_mode_sidebar(selection_state, gui_state, mesh_manager, simulation, physics_renderers);
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
