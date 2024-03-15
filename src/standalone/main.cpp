#include <Eigen/Dense>
#include "mandos.hpp"
#include "viewmandos.hpp"
#include "../mesh.hpp"
#include "../async_simulation_loop.hpp"


int main(void) {
    const Scalar minAngle = 0;
    const Scalar maxAngle = 2*M_PI;
    const int samples = 1000;
    const Vec3 axis = Vec3::Random().normalized();

    for (int i = 0; i < samples; i++) {
        const Scalar angle = (maxAngle - minAngle) * static_cast<Scalar>(i) / static_cast<Scalar>(samples);
        Eigen::Matrix<Scalar, 3, 9> dvecR_dphi = dvecR_dtheta_global(angle * axis);
        Eigen::JacobiSVD<Eigen::Matrix<Scalar,3, 9>> global_svd(dvecR_dphi);
        std::cout << angle << " ";
        std ::cout << global_svd.singularValues().transpose() << std ::endl;
    }

    // // Geometry loading
    // RenderMesh renderMesh = RenderMesh("resources/obj/triceratops3.obj");
    // SimulationMesh simMesh = SimulationMesh(renderMesh);
    // TetrahedronMesh tmesh = TetrahedronMesh(simMesh);

    // // Simulation description
    // Simulation simulation;
    // const Scalar mass = 1000;
    // const Scalar poisson_ratio = 0.2;
    // const Scalar young_modulus = 5e3;

    // FEMHandle fem = FEMHandle(simulation, tmesh.vertices, tmesh.indices,
    //                           mass, poisson_ratio, young_modulus)
    //     .add_gravity(-1.0)
    //     ;
    // ParticleHandle p1 = ParticleHandle(simulation, mass).freeze();
    // ParticleHandle p2 = get_particle_handle(simulation, 0);
    // Vec3 pos = p2.get_position(simulation.initial_state);
    // p1.set_initial_position(pos + Vec3(0,1,0));

    // join_particles_with_spring(simulation, p1, p2, 1000.0, 0.0);
    // PhysicsState state = simulation.initial_state;

    // // Start the simulation thread
    // simulation_async_loop(simulation);

    // // Render
    // MandosViewer viewer;
    // MeshGPU meshGPU = MeshGPU(renderMesh);

    // // Render loop
    // bool simulation_paused = false;
    // while (not viewer.window_should_close()) {

    //     // Interaction with the simulaiton
    //     if (viewer.is_key_pressed(Key_H)) {
    //         simulation_paused = not simulation_paused;
    //     }

    //     if (viewer.is_key_pressed(Key_G)) {
    //         simulation_async_loop_request_iteration();
    //         state = get_current_physics_state();
    //     }

    //     if (viewer.is_key_pressed(Key_R)) {
    //         state = simulation.initial_state;
    //         set_current_physics_state(state);
    //     }

    //     if (not simulation_paused) {
    //         EnergyAndDerivatives f = EnergyAndDerivatives(0);
    //         simulation_async_loop_request_iteration();
    //         state = get_current_physics_state();
    //     }

    //     // Rendering

    //     viewer.begin_drawing();
    //     viewer.draw_simulation_state(simulation, state);
    //     viewer.draw_FEM(fem, state, meshGPU, renderMesh, simMesh);
    //     viewer.draw_particle(p1, state);
    //     viewer.draw_particle(p2, state);
    //     viewer.end_drawing();
    // }
}
