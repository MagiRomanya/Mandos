#include <Eigen/Dense>
#include "mandos.hpp"
#include "viewmandos.hpp"
#include "../rod_segment.hpp"

int main(void) {
    // Simulation description
    Simulation simulation;

    PlaneCollider plane = {
    .center = Vec3(0.0, -1.0, 0.0),
    .normal = Vec3(0.0, 1.0, 0.0),
    };
    simulation.colliders.plane_colliders.push_back(plane);

    simulation.TimeStep = 0.05;
    const Scalar mass = 10;
    const Scalar gravity = -9.8;

    const unsigned int nRigidBodies = 30;
    const Scalar L0 = 0.5;
    const Scalar midpoint = L0 * nRigidBodies * 0.5;
    const Vec3 fixer_pos0 = Vec3(0, 4, 0);
    const Vec3 fixer_pos1 = Vec3(0, 1, -midpoint);
    const Vec3 fixer_pos2 = Vec3(0, 1, midpoint - L0);
    ParticleHandle p_fix0 = ParticleHandle(simulation,1)
        .freeze()
        .set_initial_position(fixer_pos0);
    ParticleHandle p_fix1 = ParticleHandle(simulation,1)
        .freeze()
        .set_initial_position(fixer_pos1);

    ParticleHandle p_fix2 = ParticleHandle(simulation,1)
        .freeze()
        .set_initial_position(fixer_pos2);

    const RodSegmentParameters parameters = {
    .Ks = 5000.0,
    .L0 = L0,
    .translational_damping = 5000.0,
    .rotational_damping = 0.0,
    .constraint_stiffness = 6000.0,
    .intrinsic_darboux = Vec3::Zero(),
    // .intrinsic_darboux = Vec3(0.0, 0.0, 0.1),
    .stiffness_tensor = 5000.0 * Vec3::Ones(),
    // .stiffness_tensor = Vec3(500000.0, 500000.0, 500000.0),
    };

    // generate_rod(simulation, nRigidBodies-1, 10*nRigidBodies, 15.0, Vec3(0,0,0), Vec3(0,0,1), parameters);
    RodHandle rod = RodHandle(simulation, nRigidBodies, 15.0, 10*nRigidBodies, parameters)
        .set_initial_origin_position(Vec3(-7.5 ,0.0 ,0.0 ))
        .set_initial_rod_direction(Vec3(1.0, 1.0, 0.0))
        .add_gravity(gravity)
        ;


    RigidBody rb = simulation.simulables.rigid_bodies[0];
    RigidBody rb2 = simulation.simulables.rigid_bodies[nRigidBodies-1];

    // for (unsigned int i = 0; i < 6; i++) {
    //     simulation.frozen_dof.push_back(rb.index+i);
    //     // simulation.frozen_dof.push_back(rb2.index+i);
    // }

    SpringParameters spring_param0 = SpringParameters(500.0,
                                                     (rb.get_COM_position(simulation.initial_state.x) - fixer_pos0).norm(),
                                                     10.0);
    SpringParameters spring_param1 = SpringParameters(500.0,
                                                     (rb.get_COM_position(simulation.initial_state.x) - fixer_pos0).norm(),
                                                     10.0);
    SpringParameters spring_param2 = SpringParameters(500.0,
                                                     (rb2.get_COM_position(simulation.initial_state.x) - fixer_pos0).norm(),
                                                     10.0);

    // simulation.energies.particle_springs.emplace_back(Particle(rb.mass, rb.index), p_fix0.particle, spring_param1);
    // simulation.energies.particle_springs.emplace_back(Particle(rb2.mass, rb2.index), p_fix0.particle, spring_param2);

    // simulation.energies.particle_springs.emplace_back(Particle(rb.mass, rb.index), p_fix1.particle, spring_param1);
    // simulation.energies.particle_springs.emplace_back(Particle(rb2.mass, rb2.index), p_fix2.particle, spring_param2);


    PhysicsState state = simulation.initial_state;
    // Render loop
    MandosViewer viewer;
    bool simulation_paused = true;
    Scalar time = 0.0;
    while (not viewer.window_should_close()) {

        // Interaction with the simulaiton
        if (viewer.is_key_pressed(Key_H)) {
            simulation_paused = not simulation_paused;
        }

        if (viewer.is_key_pressed(Key_G)) {
            simulation_step(simulation, state);
            time += simulation.TimeStep;
        }

        if (viewer.is_key_pressed(Key_R)) {
            state = simulation.initial_state;
            time = 0.0;
        }

        if (not simulation_paused) {
            EnergyAndDerivatives f = EnergyAndDerivatives(0);
            simulation_step(simulation, state, f);
            time += simulation.TimeStep;
            DEBUG_LOG(f.energy);
        }

        // Vec3 rotation_vector = Vec3(0, 0, 1) * std::fmod( 0.1 * time, 2 * M_PI);
        // state.x.segment<3>(rb.index + 3) = rotation_vector;

        // Vec3 position = simulation.initial_state.x.segment<3>(rb.index) + Vec3(std::sin(2*time),std::cos(2*time) - 1.0, 0);
        // // Vec3 position = simulation.initial_state.x.segment<3>(rb.index) + Vec3(0, 0, 0.1*std::sin(2 * time));
        // state.x.segment<3>(rb.index) =  position;

        // Rendering
        viewer.begin_drawing();
        viewer.draw_simulation_state(simulation, state);

        // viewer.draw_particles(simulation, state);
        viewer.draw_rigid_bodies(simulation, state);
        viewer.draw_springs(simulation, state);
        // viewer.draw_rods(simulation, state);
        // viewer.draw_colliders(simulation);

        viewer.end_drawing();
    }
}
