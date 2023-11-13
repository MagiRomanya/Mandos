#include <algorithm>
#include <imgui.h>
#include <string>
#include <unordered_map>
#include <raylib.h>

#include "gui/edit_mode.hpp"

std::unordered_map<std::string, Mesh> meshes;

void edit_mode_render_meshes() {
    for (const auto& pair: meshes) {
        DrawModel(LoadModelFromMesh(pair.second), Vector3{0}, 1, BLUE);
    }
}

void edit_mode_sidebar() {
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
}
