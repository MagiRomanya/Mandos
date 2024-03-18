#include <Eigen/Dense>
#include "mandos.hpp"
#include "viewmandos.hpp"
#include "../mesh.hpp"
#include "../async_simulation_loop.hpp"


inline Eigen::Matrix<Scalar,3,9> dvecR_dtheta_local_matrix(const Mat3& R) {
    return vectorized_levi_civita() * block_matrix(R);
}

inline Eigen::Matrix<Scalar,3,9> dvecR_dtheta_local_finite(const Mat3& R) {
    const Scalar dx = 1e-6;
    Eigen::Matrix<Scalar,3,9> dvecR_dthtea;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 theta = Vec3::Zero();
        theta(i) = dx;
        const Mat3 newR = compute_rotation_matrix_rodrigues(theta) * R;
        dvecR_dthtea.row(i) = (vectorize_matrix(newR) - vectorize_matrix(R)) / dx;
    }
    return dvecR_dthtea;
}

Vec3 compute_darboux_vector(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Mat3 dR_dx = (R2 - R1) / L0;
    const Mat3 R = (R2 + R1) / 2;
    const Mat3 skew_u = dR_dx * R.transpose();
    const Vec3 u = 0.5 * vectorized_levi_civita() * vectorize_matrix<3>(skew_u);
    // const Vec3 u_2 = Vec3(-skew_u(1,2), skew_u(0,2), -skew_u(0,1));;
    // DEBUG_LOG(u.transpose());
    // DEBUG_LOG(u_2.transpose());
    return u;
}

Mat3 compute_darboux_vector_theta1_derivative(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Mat3 dR_dx = (R2 - R1) / L0;
    const Mat3 R = (R2 + R1) / 2;

    Eigen::Matrix<Scalar,3,9> dvecR_dtheta = 0.5 * dvecR_dtheta_local_matrix(R2);
    const Eigen::Matrix<Scalar,3,9> dvecRx_dtheta = dvecR_dtheta_local_matrix(R2) / L0;

    const Eigen::Matrix<Scalar,9,3> dvec_skewU_dtheta =
        block_matrix<3,3>(R) * dvecRx_dtheta.transpose()                                           // d(dR_dx)/dtheta RT
        + transpose_vectorized_matrix_N<9,3>(block_matrix<3,3>(dR_dx) * dvecR_dtheta.transpose())   // dR/dx (dRT/dtheta)
        ;
    Mat3 du_dtheta = 0.5 * vectorized_levi_civita() * dvec_skewU_dtheta;
    return du_dtheta;
}

Mat3 compute_darboux_vector_theta2_derivative(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Mat3 dR_dx = (R2 - R1) / L0;
    const Mat3 R = (R2 + R1) / 2;

    Eigen::Matrix<Scalar,3,9> dvecR_dtheta = 0.5 * dvecR_dtheta_local_matrix(R1);
    const Eigen::Matrix<Scalar,3,9> dvecRx_dtheta = - dvecR_dtheta_local_matrix(R1) / L0;

    const Eigen::Matrix<Scalar,9,3> dvec_skewU_dtheta =
        block_matrix<3,3>(R) * dvecRx_dtheta.transpose()                                           // d(dR_dx)/dtheta RT
        + transpose_vectorized_matrix_N<9,3>(block_matrix<3,3>(dR_dx) * dvecR_dtheta.transpose())   // dR/dx (dRT/dtheta)
        ;
    Mat3 du_dtheta = 0.5 * vectorized_levi_civita() * dvec_skewU_dtheta;
    return du_dtheta;
}


Mat3 compute_darboux_vector_local_finite_derivative(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Vec3 u0 = compute_darboux_vector(L0, R1, R2);
    const Scalar dx = 1e-6;
    Mat3 du_dtheta;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 theta = Vec3::Zero();
        theta(i) = dx;
        // const Mat3 newR2 = compute_rotation_matrix_rodrigues(theta) * R2;
        // const Vec3 newU = compute_darboux_vector(L0, R1, newR2);
        const Mat3 newR1 = compute_rotation_matrix_rodrigues(theta) * R1;
        const Vec3 newU = compute_darboux_vector(L0, newR1, R2);

        du_dtheta.col(i) = (newU - u0) / dx;
    }
    return du_dtheta;
}

int main(void) {
    const Mat3 R1 = compute_rotation_matrix_rodrigues(Vec3(1,1,1).normalized());
    const Mat3 R2 = compute_rotation_matrix_rodrigues(Vec3(1,-1,-1).normalized());
    const Scalar L0 = 1.0;
    // const Vec3 u = compute_darboux_vector(L0, R1, R2);

    const Mat3 du_dtheta_finite = compute_darboux_vector_local_finite_derivative(L0, R1, R2);
    const Mat3 du_dtheta = compute_darboux_vector_theta1_derivative(L0, R1, R2);
    const Mat3 du_dtheta2 = compute_darboux_vector_theta2_derivative(L0, R1, R2);
    std ::cout << "du_dtheta_finite"
               << "\n" << du_dtheta_finite << std ::endl;
    std ::cout << "du_dtheta"
               << "\n" << du_dtheta << std ::endl;
    std ::cout << "du_dtheta2"
               << "\n" << du_dtheta2 << std ::endl;

    // const Scalar minAngle = 0;
    // const Scalar maxAngle = M_PI;
    // const int samples = 1000;
    // const Scalar magnitude = 1;
    // const Vec3 axis = Vec3::Random().normalized();

    // for (int i = 0; i < samples; i++) {
    //     const Scalar angle = (maxAngle - minAngle) * static_cast<Scalar>(i) / static_cast<Scalar>(samples);
    //     // Eigen::Matrix<Scalar, 3, 9> dvecR_dphi = dvecR_dtheta_global(angle * axis);
    //     // Eigen::JacobiSVD<Eigen::Matrix<Scalar,3, 9>> global_svd(dvecR_dphi);

    //     Eigen::EigenSolver<Mat3> solver(compute_global_to_local_axis_angle_jacobian(axis * angle));
    //     std::cout << angle << " ";
    //     std ::cout << solver.eigenvalues().real().transpose() << " ";
    //     std ::cout << solver.eigenvalues().imag().transpose() << std ::endl;
    //     // std ::cout << global_svd.singularValues().transpose() << std ::endl;
    // }

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
