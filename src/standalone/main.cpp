#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

#include "linear_algebra.hpp"
#include "mandos.hpp"
#include "mesh.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "clock.hpp"
#include "viewmandos.hpp"


int main(void) {
    Simulation simulation;
    const Scalar MASS = 1.0;
    const Scalar GRAVITY = -1.0;

    RenderMesh render_mesh = RenderMesh("resources/obj/cone.obj");
    SimulationMesh RB_mesh = SimulationMesh(render_mesh);

    const Mat3 inertia_tensor = compute_initial_inertia_tensor_PARTICLES(MASS, RB_mesh.vertices);
    const Vec3 center_of_mass = compute_COM_position_PARTICLES(RB_mesh.vertices);

    recenter_mesh(RB_mesh, center_of_mass);


    const RigidBodyHandle rb1 = RigidBodyHandle(simulation, MASS, inertia_tensor)
        .add_gravity(GRAVITY)
        .set_COM_initial_position(Vec3(0,0,0));

    const ParticleHandle p1 = ParticleHandle(simulation, MASS/10.0)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(1,2,0));

    const ParticleHandle p2 = ParticleHandle(simulation, MASS/10.0)
        .freeze()
        .set_initial_position(Vec3(0,0,1));

    const ParticleHandle p3 = ParticleHandle(simulation, MASS/10.0)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(0,1,-1));

    const ParticleHandle p4 = ParticleHandle(simulation, MASS/10.0)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(0,1,1));

    const ParticleHandle p5 = ParticleHandle(simulation, MASS/10.0)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(1,1,1));

    const ParticleHandle p6 = ParticleHandle(simulation, MASS/10.0)
        .add_gravity(GRAVITY)
        .set_initial_position(Vec3(2,1,1));

    std::vector<ParticleRigidBodyCopuling> copulings_vec;
    ParticleRigidBodyCopuling copuling = ParticleRigidBodyCopuling(rb1.rb, p1.particle, p1.particle.get_position(simulation.initial_state.x));
    copulings_vec.push_back(copuling);
    ParticleRigidBodyCopuling copuling2 = ParticleRigidBodyCopuling(rb1.rb, p3.particle, p3.particle.get_position(simulation.initial_state.x));
    copulings_vec.push_back(copuling2);
    ParticleRigidBodyCopuling copuling3 = ParticleRigidBodyCopuling(rb1.rb, p4.particle, p4.particle.get_position(simulation.initial_state.x));
    copulings_vec.push_back(copuling3);

    simulation.copulings = Copulings(copulings_vec);

    // std::cout << rb1.rb.get_COM_position(simulation.initial_state.x).transpose() << std::endl;
    // std::cout << rb1.rb.get_COM_position(simulation.initial_state.x_old).transpose() << std::endl;
    // std::cout << simulation.initial_state.x.transpose() << std::endl;

    join_particles_with_spring(simulation, p1, p2, 10.0, 1);
    join_particles_with_spring(simulation, p4, p5, 10.0, 1);
    join_particles_with_spring(simulation, p5, p6, 10.0, 1);

    //--------------------------------------------------------------------------------------

    bool simulation_paused = true;
    PhysicsState state = simulation.initial_state;

    // Render Initialization
    //--------------------------------------------------------------------------------------
    MandosViewer viewer = MandosViewer();
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
            // std::cout << "Energy " << f.energy << std::endl;
            std::cout << "Simulation time " << simulation_time << " (ms)" << std::endl;

        }
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        viewer.begin_drawing();
        {
            viewer.begin_3D_mode();
            {
                viewer.draw_rigid_body(rb1, state, meshGPU);
                viewer.draw_particle(p1, state);
                viewer.draw_particle(p2, state);
                viewer.draw_particle(p3, state);
                viewer.draw_particle(p4, state);
                viewer.draw_particle(p5, state);
                viewer.draw_particle(p6, state);
            }
            viewer.end_3D_mode();
        }
        viewer.end_drawing();
        //----------------------------------------------------------------------------------
    }

    return 0;
}
