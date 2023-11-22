#include <iostream>
#include <vector>
#include <raylib.h>
#include <rlgl.h>


#include "fem_unit.hpp"
#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "render/simulation_visualization.hpp"
#include "simulation.hpp"

int main(void) {
    // Create tetrahedron

    Simulation simulation;
    const Scalar particle_mass = 1;

    // Add 3 dimensions per particle to the degrees of freedom
    const unsigned int index = simulation.initial_state.x.size();
    simulation.initial_state.add_size(4*3);

    // Create particle simulables
    simulation.simulables.particles.emplace_back(particle_mass, index+3*0);
    simulation.simulables.particles.emplace_back(particle_mass, index+3*1);
    simulation.simulables.particles.emplace_back(particle_mass, index+3*2);
    simulation.simulables.particles.emplace_back(particle_mass, index+3*3);

    // Initial conditions
    simulation.initial_state.x.segment<3>(3*0+index) = Vec3(0,0,0);
    simulation.initial_state.x.segment<3>(3*1+index) = Vec3(1,0,0);
    simulation.initial_state.x.segment<3>(3*2+index) = Vec3(0,1,0);
    simulation.initial_state.x.segment<3>(3*3+index) = Vec3(0,0,1);
    simulation.initial_state.v.setZero();
    // const Scalar vel = 1;
    // simulation.initial_state.v.segment<3>(3*2+index) = Vec3(0, vel/2,0);
    // simulation.initial_state.v.segment<3>(3*0+index) = Vec3(0,-vel/2,0);

    const Scalar young_modulus = 50e3;
    const Scalar poisson_ratio = 0.2;
    const Scalar mu = young_modulus / (2*(1 + poisson_ratio));
    const Scalar lambda = young_modulus * poisson_ratio / ((1+poisson_ratio) * (1-2*poisson_ratio));
    Eigen::Matrix<Scalar,4,3> ds_dx =  compute_shape_function_derivative(Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1));
    FEM_ElementParameters FEM_paramters(mu, lambda, ds_dx);
    simulation.energies.fem_elements_3d.emplace_back(simulation.simulables.particles[0],
                                                     simulation.simulables.particles[1],
                                                     simulation.simulables.particles[2],
                                                     simulation.simulables.particles[3],
                                                     FEM_paramters);

    // Displace original position
    simulation.initial_state.x.segment<3>(3*2+index) = Vec3(0,1.5,0);

    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos simulator");
    Camera3D camera = create_camera();

    SetTargetFPS(200);
    //--------------------------------------------------------------------------------------

    PhysicsState state = simulation.initial_state;
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera, CAMERA_ORBITAL);
        /// Keyboard controls
        if (IsKeyPressed(KEY_Q)) break;

        if (true or IsKeyPressed(KEY_G)) {
            EnergyAndDerivatives f(0);
            simulation_step(simulation, state, f);
            std::cout << "Energy " << f.energy << std::endl;
        }
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);
            {
                DrawGrid(30, 1.0f);
                for (unsigned int i = 0; i < simulation.energies.fem_elements_3d.size(); i++) {
                    const FEM_Element3D& e = simulation.energies.fem_elements_3d[i];
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
                Color colors[] = {RED, BLUE, GREEN, MAGENTA};
                for (unsigned int i = 0; i < simulation.simulables.particles.size(); i++) {
                    const Particle& p = simulation.simulables.particles[i];
                    const Vec3 x = p.get_position(state);
                    DrawSphere(Vector3{x.x(), x.y(), x.z()}, 0.05, colors[i]);
                }
            }
            EndMode3D();

            DrawFPS(10, 10);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }
    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
