#include <Eigen/Dense>
#include "viewmandos.hpp"
#include "../rod_segment.hpp"
#include "../mesh.hpp"
#include "../async_simulation_loop.hpp"

int main(void) {
    // Simulation description
    Simulation simulation;
    const Scalar mass = 1;

    const unsigned int nRigidBodies = 10;
    const Scalar L0 = 1.0;
    for (unsigned int i = 0; i < nRigidBodies; i++) {
        // Create rigid_body
        unsigned int index = simulation.initial_state.get_nDoF();

        RigidBody rb = RigidBody(index, mass, Mat3::Identity());
        add_rigid_body_to_simulation(simulation, rb);

        // Initial state
        simulation.initial_state.x.segment<6>(index) = Vec6(0.0, 0.0, L0*i, 0.0, 0.0, 0.0);
        simulation.initial_state.v.segment<6>(index) = Vec6(0.0, 0.0,  0.0, 0.0, 0.0, 0.0); // At rest
    }

    // Create the rod segments
    const RodSegmentParameters parameters = {
    .Ks = 10.0,
    .L0 = L0,
    .translational_damping = 0.0,
    .rotational_damping = 0.0,
    .constraint_stiffness = 10.0,
    .intrinsic_darboux = Vec3::Zero(),
    .stiffness_tensor = 10.0 * Vec3::Ones(),
    };

    for (unsigned int i = 0; i < nRigidBodies - 1; i++) {
        const RigidBody rbA = simulation.simulables.rigid_bodies[i];
        const RigidBody rbB = simulation.simulables.rigid_bodies[i+1];
        RodSegment segment = RodSegment(rbA, rbB, parameters);
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
            simulation_async_loop_request_iteration();
            state = get_current_physics_state();
        }

        if (viewer.is_key_pressed(Key_R)) {
            state = simulation.initial_state;
            set_current_physics_state(state);
        }

        if (not simulation_paused) {
            EnergyAndDerivatives f = EnergyAndDerivatives(0);
            simulation_async_loop_request_iteration();
            state = get_current_physics_state();
        }

        // Rendering

        viewer.begin_drawing();
        viewer.draw_simulation_state(simulation, state);

        viewer.draw_rigid_bodies(simulation, state);

        viewer.end_drawing();
    }
}
