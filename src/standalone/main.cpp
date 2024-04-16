#include <Eigen/Dense>
#include "mandos.hpp"
#include "viewmandos.hpp"
#include "../rod_segment.hpp"

int main(void) {
    // Simulation description
    Simulation simulation;

    SimulationMesh collider_mesh = SimulationMesh("resources/obj/templateM.obj");
    // SimulationMesh collider_mesh = SimulationMesh("resources/obj/monke.obj");
    SDFColliderHandle sdf_collider = SDFColliderHandle(simulation, collider_mesh);
    PlaneColliderHandle plane_collider = PlaneColliderHandle(simulation)
        .set_direction(Vec3(0.0, 1.0, 0.0))
        .set_origin_position(Vec3(0.0, -30.0, 0.0))
        ;

    simulation.TimeStep = 0.05;
    const Scalar mass = 10;
    const Scalar gravity = -9.8 * mass;

    const RodSegmentParameters parameters = {
    .Ks = 500000.0,
    .L0 = 0.5,
    .translational_damping = 5000.0,
    .rotational_damping = 0.0,
    .constraint_stiffness = 20000.0,
    .intrinsic_darboux = Vec3::Zero(),
    // .stiffness_tensor = 500000000.0 * Vec3::Ones(),
    .stiffness_tensor = 5000.0 * Vec3::Ones(),
    // .stiffness_tensor = Vec3(500000.0, 500000.0, 500000.0),
    };


    const Scalar height = 20.0;
    const Scalar length = 40.0;
    // const Scalar length = 10.0;
    const unsigned int segments = 100;
    RodHandle rod1 = RodHandle(simulation, segments, length*Vec3(1.0, 0.0, 0.0).normalized(), mass * segments, parameters)
        .set_initial_rod_position(Vec3(-length * 0.5, height, 0.0 ))
        .add_gravity(gravity)
        ;

    RodHandle rod2 = RodHandle(simulation, segments, length*Vec3(0.0, 0.0, 1.0).normalized(), mass*segments, parameters)
        .set_initial_rod_position(Vec3(0, height , -length * 0.5))
        .add_gravity(gravity)
        ;

    RodHandle rod3 = RodHandle(simulation, segments, length*Vec3(1.0, 0.0, 1.0).normalized(), mass*segments, parameters)
        .set_initial_rod_position(Vec3(1, 0, 1).normalized() * (-length * 0.5) + Vec3(0,height,0))
        .add_gravity(gravity)
        ;

    RodHandle rod4 = RodHandle(simulation, segments, length*Vec3(1.0, 0.0, -1.0).normalized(), mass*segments, parameters)
        .set_initial_rod_position(Vec3(1, 0, -1).normalized() * (-length * 0.5) + Vec3(0,height,0))
        .add_gravity(gravity)
        ;



    // std::vector<Scalar> spring_curve;
    // LoadCurveTinyObj("resources/obj/spring-curve.obj", spring_curve);
    // RodHandle rod = RodHandle(simulation, spring_curve, 100, parameters)
    //     .set_initial_rod_position(Vec3(0, 20.0, -4.0))
    //     .add_gravity(gravity)
    //     ;


    // RigidBody rb = simulation.simulables.rigid_bodies[0];
    // RigidBody rb2 = simulation.simulables.rigid_bodies[rod.get_n_rigid_bodies()-1];

    // for (unsigned int i = 0; i < 6; i++) {
    //     // simulation.frozen_dof.push_back(rb.index+i);
    //     simulation.frozen_dof.push_back(rb2.index+i);
    // }


    PhysicsState state = simulation.initial_state;
    // Render loop
    MandosViewer viewer;
    RenderMesh rodMesh = RenderMesh("resources/obj/skinning-cylinder.obj");
    SkinnedRodGPU rodGPU = SkinnedRodGPU(rodMesh, rod1.get_n_rigid_bodies() - 1);
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
        // viewer.draw_rigid_bodies(simulation, state);
        viewer.draw_springs(simulation, state);
        viewer.draw_rod(rod1, state, rodGPU);
        viewer.draw_rod(rod2, state, rodGPU);
        viewer.draw_rod(rod3, state, rodGPU);
        viewer.draw_rod(rod4, state, rodGPU);
        // viewer.draw_rods(simulation, state);

        viewer.end_drawing();
    }

    Eigen::Matrix<Scalar,2,Eigen::Dynamic> boneWeights;
    Eigen::Matrix<unsigned int, 2, Eigen::Dynamic> boneIDs;
}
