#include <Eigen/Dense>
#include "mandos.hpp"
#include "viewmandos.hpp"
#include "../rod_segment.hpp"
#include "../mesh.hpp"
#include "../async_simulation_loop.hpp"

int main(void) {
    // Simulation description
    Simulation simulation;
    simulation.TimeStep = 0.05;
    const Scalar mass = 10;
    const Scalar gravity = -1;

    const unsigned int nRigidBodies = 30;
    const Scalar L0 = 0.5;
    const Scalar midpoint = L0 * nRigidBodies * 0.5;
    const Vec3 fixer_pos1 = Vec3(0, 1, -midpoint);
    const Vec3 fixer_pos2 = Vec3(0, 1, midpoint);
    ParticleHandle p_fix1 = ParticleHandle(simulation,1)
        .freeze()
        .set_initial_position(fixer_pos1);

    ParticleHandle p_fix2 = ParticleHandle(simulation,1)
        .freeze()
        .set_initial_position(fixer_pos2);

    for (unsigned int i = 0; i < nRigidBodies; i++) {
        // Create rigid_body
        unsigned int index = simulation.initial_state.get_nDoF();
        simulation.initial_state.add_size(6);
        RigidBody rb = RigidBody(index, mass, mass * Mat3::Identity());
        add_rigid_body_to_simulation(simulation, rb);

        // Initial state
        simulation.initial_state.x.segment<6>(index) = Vec6(0.0, 0.0, L0*i - midpoint, 0.0, 0.0, 0.0);
        simulation.initial_state.v.segment<6>(index) = Vec6(0.0, 0.0,  0.0, 0.0, 0.0, 0.0); // At rest

        // Add gravity
        unsigned int y_index = index + 1;
        simulation.energies.gravities.emplace_back(y_index, GravityParameters(gravity));
    }

    // Join the first chain element to a fixed particle with a spring
    RigidBody rb = simulation.simulables.rigid_bodies[0];
    RigidBody rb2 = simulation.simulables.rigid_bodies[29];
    SpringParameters spring_param = SpringParameters(50.0,
                                                     (rb.get_COM_position(simulation.initial_state.x) - fixer_pos1).norm(),
                                                     0);
    SpringParameters spring_param2 = SpringParameters(50.0,
                                                     (rb.get_COM_position(simulation.initial_state.x) - fixer_pos1).norm(),
                                                     0);
    simulation.energies.particle_springs.emplace_back(Particle(rb.mass, rb.index), p_fix1.particle, spring_param);
    simulation.energies.particle_springs.emplace_back(Particle(rb2.mass, rb2.index), p_fix2.particle, spring_param2);

    // Create the rod segments
    const RodSegmentParameters parameters = {
    .Ks = 200.0,
    .L0 = L0,
    .translational_damping = 10.0,
    .rotational_damping = 0.0,
    .constraint_stiffness = 500.0,
    .intrinsic_darboux = Vec3::Zero(),
    .stiffness_tensor = 100.0 * Vec3::Ones(),
    };

    for (unsigned int i = 0; i < nRigidBodies - 1; i++) {
        const RigidBody rbA = simulation.simulables.rigid_bodies[i];
        const RigidBody rbB = simulation.simulables.rigid_bodies[i+1];
        RodSegment segment = RodSegment(rbB, rbA, parameters);

        simulation.energies.rod_segments.push_back(segment);
    }


    PhysicsState state = simulation.initial_state;
    // Render loop
    MandosViewer viewer;
    bool simulation_paused = true;
    while (not viewer.window_should_close()) {

        // Interaction with the simulaiton
        if (viewer.is_key_pressed(Key_H)) {
            simulation_paused = not simulation_paused;
        }

        if (viewer.is_key_pressed(Key_G)) {
            simulation_step(simulation, state);
        }

        if (viewer.is_key_pressed(Key_R)) {
            state = simulation.initial_state;
        }

        if (not simulation_paused) {
            EnergyAndDerivatives f = EnergyAndDerivatives(0);
            simulation_step(simulation, state, f);
            DEBUG_LOG(f.energy);
        }

        // Rendering
        viewer.begin_drawing();
        viewer.draw_simulation_state(simulation, state);

        viewer.draw_particles(simulation, state);
        viewer.draw_rigid_bodies(simulation, state);
        viewer.draw_springs(simulation, state);

        viewer.end_drawing();
    }
}
