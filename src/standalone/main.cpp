#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>
#include <rlImGui.h>
#include <imgui.h>

#include "linear_algebra.hpp"
#include "physics_render.hpp"
#include "physics_state.hpp"
#include "render/simulation_visualization.hpp"
#include "simulation.hpp"
#include "simulable_generator.hpp"
#include "integrators.hpp"
#include "utility_functions.hpp"

int main(int argc, char *argv[]) {
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "raylib [core] example - 3d camera mode");
    rlImGuiSetup(false);
    // Define the camera to look into our 3d world
    rlDisableBackfaceCulling();
    Camera3D camera = { 0 };
    camera.position = Vector3{ 0.0f, 4.0f, 10.0f };  // Camera position
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type

    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------

    Mesh cloth_mesh = GenMeshPlane(10.0, 10.0, 20, 20);
    const unsigned int nDoF = cloth_mesh.vertexCount * 3;
    const unsigned int n_indices = cloth_mesh.triangleCount * 3;
    std::vector<Scalar> vertices(cloth_mesh.vertices, cloth_mesh.vertices + nDoF);
    std::vector<unsigned int> indices(cloth_mesh.indices, cloth_mesh.indices + n_indices);
    std::cout << "Number of degrees of freedom: " << nDoF << std::endl;

    Simulation simulation;
    SimulableBounds cloth_bounds = generate_mass_spring(simulation, vertices, indices, 0.1, 100.0, 1.0);
    MassSpringRenderer cloth_renderer = MassSpringRenderer(cloth_mesh, cloth_bounds);
    PhysicsRenderers renderers;
    renderers.mass_spring.push_back(cloth_renderer);

    // froze degrees of freedom
    simulation.frozen_dof.push_back(0);
    simulation.frozen_dof.push_back(1);
    simulation.frozen_dof.push_back(2);

    {
        // Test volume equation
        Mesh mesh = LoadMeshTinyOBJ("img/obj/sphere.obj");
        const unsigned int n_vertices = mesh.vertexCount * 3;
        const unsigned int n_indices = mesh.triangleCount * 3;
        const std::vector<Scalar> vertices(mesh.vertices, mesh.vertices + n_vertices);
        const std::vector<unsigned int> indices(mesh.indices, mesh.indices + n_indices);
        const Scalar volume = compute_mesh_volume(indices, vertices);
        std::cout << "Sphere volume " << volume << std::endl;
    }

    // simulation_visualization_loop(simulation, renderers);

    // De-Initialization
    //--------------------------------------------------------------------------------------
    rlImGuiShutdown();
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
