#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

#include "linear_algebra.hpp"
#include "mandos.hpp"
#include "mesh.hpp"
#include "rigid_body.hpp"
#include "clock.hpp"
#include "utility_functions.hpp"
#include "viewmandos.hpp"


int main(void) {
    Simulation simulation;
    const Scalar MASS = 10.0;
    const Scalar GRAVITY = -1.0;

    RenderMesh render_mesh = RenderMesh("resources/obj/bunny.obj");
    SimulationMesh sim_mesh = SimulationMesh("resources/obj/bunny.obj");

    const Mat3 inertia_tensor = compute_initial_inertia_tensor_PARTICLES(MASS, sim_mesh.vertices);
    const Vec3 center_of_mass = compute_COM_position_PARTICLES(sim_mesh.vertices);
    recenter_mesh(sim_mesh, center_of_mass);

    // Compute the tetrahedron mesh
    std::vector<unsigned int> tet_indices;
    std::vector<float> tet_vertices;
    tetgen_compute_tetrahedrons(sim_mesh.indices, sim_mesh.vertices, tet_indices, tet_vertices);

    const Scalar young_modulus = 5e2;
    const Scalar poisson_ratio = 0.2;
    const FEMHandle fem1 = FEMHandle(simulation, tet_vertices, tet_indices, MASS, poisson_ratio, young_modulus)
        .add_gravity(GRAVITY)
        // .freeze_particles({0,1,4,5})
        .freeze_particles({0})
        ;
    //--------------------------------------------------------------------------------------

    bool simulation_paused = true;
    PhysicsState state = simulation.initial_state;

    // Render Initialization
    //--------------------------------------------------------------------------------------
    MandosViewer viewer = MandosViewer();
    render_mesh.updateFromSimulationMesh(sim_mesh);
    MeshGPU meshGPU = MeshGPU(render_mesh);

    while (!viewer.window_should_close())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        viewer.update_camera();
        /// Keyboard controls
        if (viewer.is_key_pressed(Key_Q)) break;

        if (viewer.is_key_pressed(Key_H)) simulation_paused = !simulation_paused;
        if (viewer.is_key_pressed(Key_R)) state = simulation.initial_state;
        static double simulation_time = 0;
        if (!simulation_paused or viewer.is_key_pressed(Key_G)) {
            Clock clock = Clock(simulation_time);
            EnergyAndDerivatives f(0);
            simulation_step(simulation, state, f);
            std::cout << "Energy " << f.energy << std::endl;
            // std::cout << "Simulation time " << simulation_time << " (ms)" << std::endl;
        }
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        viewer.begin_drawing();
        {
            viewer.begin_3D_mode();
            {
                // viewer.draw_particles(simulation, state);
                // viewer.draw_FEM_tetrahedrons(simulation, state);
                viewer.draw_FEM(fem1, state, meshGPU, render_mesh, sim_mesh);
                // viewer.draw_mesh(Mat4::Identity(), meshGPU);
            }
            viewer.end_3D_mode();
        }
        viewer.end_drawing();
        //----------------------------------------------------------------------------------
    }

    return 0;
}
