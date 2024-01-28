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

Vec3 rotation_inertia_energy_gradient(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 Romega = compute_rotation_matrix_rodrigues(TimeStep * omega0);
    // const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0
    const Mat3 R_guess = Romega * R0;

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    const Mat3 A = (rot_inertia - rot_inertia.transpose()) / 2;
    const Scalar h2 = TimeStep * TimeStep;

    const Vec3 gradient = 2.0f * Vec3(-A(1,2), A(0,2), -A(0,1)) / h2;
    return gradient;
}

Mat3 rotation_inertia_energy_hessian(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 Romega = compute_rotation_matrix_rodrigues(TimeStep * omega0);
    // const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0
    const Mat3 R_guess = Romega * R0;

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    const Mat3 S = (rot_inertia + rot_inertia.transpose()) / 2; // Exact hessian
    // const Mat3 S = R * J_inertia_tensor * R.transpose(); // Linear approximation
    const Scalar h2 = TimeStep * TimeStep;

    const Mat3 hessian = 1.0f / h2 * (S.trace() * Mat3::Identity() - S);
    return hessian;
}

Vec3 rotation_inertia_energy_gradient2(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 Romega = compute_rotation_matrix_rodrigues(TimeStep * omega0);
    // const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0
    const Mat3 R_guess = Romega * R0;

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    DEBUG_LOG(rot_inertia.transpose());
    std ::cout << "vectorized inertia"
               << " "
               << vectorize_matrix<3>(R_guess).transpose() *
                  block_matrix<3,3>(J_inertia_tensor * R.transpose())
               << std ::endl;
    // const Mat3 A = (rot_inertia - rot_inertia.transpose()) / 2;
    const Mat3 A = (-rot_inertia.transpose()) / 2;
    const Scalar h2 = TimeStep * TimeStep;

    const Vec3 gradient = 1.0f/h2 * vectorized_levi_civita() * vectorize_matrix<3>(A);
    return gradient;
}

Mat3 rotation_inertia_finite_dgradE_dtheta0(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    Mat3 H;
    const Scalar dx = 0.0001f;
    const Vec3 grad0 = rotation_inertia_energy_gradient2(theta, theta0, omega0, TimeStep);
    for (unsigned i = 0; i < 3; i++) {
        Vec3 dtheta0 = theta0;
        dtheta0[i] += dx;
        const Vec3 grad = rotation_inertia_energy_gradient2(theta, dtheta0, omega0, TimeStep);
        H.row(i) = (grad - grad0) / dx;
    }
    return H;
}

Mat3 rotation_inertia_dgradE_dtheta0(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 Romega = compute_rotation_matrix_rodrigues(TimeStep * omega0);
    const Mat3 R_guess = Romega * R0;

    const Scalar h2 = TimeStep * TimeStep;
    const Eigen::Matrix<Scalar,3,9> vLeviCivita = vectorized_levi_civita();
    Eigen::Matrix<Scalar,9,3> dvecAdtheta0 = Eigen::Matrix<Scalar,9,3>::Zero();
    // dvecAdtheta0 += -0.5f * block_matrix<3,3>(R*J_inertia_tensor) * (vLeviCivita * block_matrix<3,3>(R_guess.transpose())).transpose();
    dvecAdtheta0 -= 0.5f * block_matrix<3,3>(J_inertia_tensor*R.transpose()) * (vLeviCivita * block_matrix<3,3>(R_guess)).transpose();
    Mat3 H = 1.0f/h2 * vLeviCivita * dvecAdtheta0;
    return H;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_finite(const Vec3& theta) {
    Eigen::Matrix<Scalar,3,9> dvecR_dtheta;
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta);
    const Eigen::Vector<Scalar,9> vecR0 = vectorize_matrix(R0);
    const Scalar dx = 0.0001f;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dtheta = Vec3::Zero();
        dtheta[i] += dx;
        const Eigen::Vector<Scalar,9> vecR = vectorize_matrix<3>(compute_rotation_matrix_rodrigues(dtheta) * R0);
        dvecR_dtheta.row(i) = (vecR - vecR0) / dx;
    }
    return dvecR_dtheta;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_finite2(const Vec3& theta) {
    Eigen::Matrix<Scalar,3,9> dvecR_dtheta;
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta);
    const Eigen::Vector<Scalar,9> vecR0 = vectorize_matrix(R0);
    const Scalar dx = 0.0001f;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dtheta = theta;
        dtheta[i] += dx;
        const Eigen::Vector<Scalar,9> vecR = vectorize_matrix<3>(compute_rotation_matrix_rodrigues(dtheta));
        dvecR_dtheta.row(i) = (vecR - vecR0) / dx;
    }
    return dvecR_dtheta;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_analytic(const Vec3& theta) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    return vectorized_levi_civita() * block_matrix(R);
}

int main(void) {
    Simulation simulation;
    const Scalar MASS = 10.0;
    const Scalar GRAVITY = -1.0;

    RenderMesh render_mesh = RenderMesh("resources/obj/bunny.obj");
    SimulationMesh sim_mesh = SimulationMesh("resources/obj/bunny.obj");

    const Mat3 inertia_tensor = compute_initial_inertia_tensor_PARTICLES(MASS, sim_mesh.vertices);
    const Vec3 center_of_mass = compute_COM_position_PARTICLES(sim_mesh.vertices);
    recenter_mesh(sim_mesh, center_of_mass);

    const Vec3 theta = Vec3(3.0,0,0);
    const Vec3 theta0 = Vec3(2,0,0);
    const Vec3 omega0 = Vec3(0,0,0);

    Eigen::Matrix<Scalar,3,9> dvecR_dtheta_F = dvecR_dtheta_finite(theta);
    Eigen::Matrix<Scalar,3,9> dvecR_dtheta_F2 = dvecR_dtheta_finite2(theta);
    Eigen::Matrix<Scalar,3,9> dvecR_dtheta_A = dvecR_dtheta_analytic(theta);
    std::cout << "dvecR_dtheta finite = "<< std::endl << dvecR_dtheta_F << std::endl;
    std::cout << "dvecR_dtheta finite 2 = "<< std::endl << dvecR_dtheta_F2 << std::endl;
    std::cout << "dvecR_dtheta analytic = "<< std::endl << dvecR_dtheta_A << std::endl;

    // const Vec3 grad1 = rotation_inertia_energy_gradient(theta, theta0, omega0, 0.1f);
    // const Vec3 grad2 = rotation_inertia_energy_gradient2(theta, theta0, omega0, 0.1f);
    // std ::cout << "grad1" << std::endl << grad1.transpose() << std ::endl;
    // std ::cout << "grad2" << std::endl << grad2.transpose() << std ::endl;

    // const Mat3 dEdtheta0 = rotation_inertia_dgradE_dtheta0(theta, theta0, omega0, 0.1f);
    // const Mat3 dEdtheta0Finite = rotation_inertia_finite_dgradE_dtheta0(theta, theta0, omega0, 0.1f);
    // const Mat3 Hess = rotation_inertia_energy_hessian(theta, theta0, omega0, 0.1f);
    // std ::cout << "dEdtheta0" << std::endl << dEdtheta0 << std ::endl;
    // std ::cout << "dEdtheta0Finite" << std::endl << dEdtheta0Finite << std ::endl;
    // std ::cout << "Hess" << std::endl << Hess << std ::endl;
    // std ::cout << "dEdtheta0 - Hess" << std::endl << dEdtheta0 - Hess << std ::endl;

    // Compute the tetrahedron mesh
    // std::vector<unsigned int> tet_indices;
    // std::vector<float> tet_vertices;
    // tetgen_compute_tetrahedrons(sim_mesh.indices, sim_mesh.vertices, tet_indices, tet_vertices);

    // const Scalar young_modulus = 5e2;
    // const Scalar poisson_ratio = 0.2;
    // const FEMHandle fem1 = FEMHandle(simulation, tet_vertices, tet_indices, MASS, poisson_ratio, young_modulus)
    //     .add_gravity(GRAVITY)
    //     // .freeze_particles({0,1,4,5})
    //     .freeze_particles({0})
    //     ;
    // //--------------------------------------------------------------------------------------

    // bool simulation_paused = true;
    // PhysicsState state = simulation.initial_state;

    // // Render Initialization
    // //--------------------------------------------------------------------------------------
    // MandosViewer viewer = MandosViewer();
    // render_mesh.updateFromSimulationMesh(sim_mesh);
    // MeshGPU meshGPU = MeshGPU(render_mesh);

    // while (!viewer.window_should_close())    // Detect window close button or ESC key
    // {
    //     // Update
    //     //----------------------------------------------------------------------------------
    //     /// Keyboard controls
    //     if (viewer.is_key_pressed(Key_Q)) break;

    //     if (viewer.is_key_pressed(Key_H)) simulation_paused = !simulation_paused;
    //     if (viewer.is_key_pressed(Key_R)) state = simulation.initial_state;
    //     static double simulation_time = 0;
    //     if (!simulation_paused or viewer.is_key_pressed(Key_G)) {
    //         Clock clock = Clock(simulation_time);
    //         EnergyAndDerivatives f(0);
    //         simulation_step(simulation, state, f);
    //         std::cout << "Energy " << f.energy << std::endl;
    //         // std::cout << "Simulation time " << simulation_time << " (ms)" << std::endl;
    //     }
    //     //----------------------------------------------------------------------------------

    //     // Draw
    //     //----------------------------------------------------------------------------------
    //     viewer.begin_drawing();
    //     {
    //         viewer.begin_3D_mode();
    //         {
    //             // viewer.draw_particles(simulation, state);
    //             // viewer.draw_FEM_tetrahedrons(simulation, state);
    //             viewer.draw_FEM(fem1, state, meshGPU, render_mesh, sim_mesh);
    //             // viewer.draw_mesh(Mat4::Identity(), meshGPU);
    //         }
    //         viewer.end_3D_mode();
    //     }
    //     viewer.end_drawing();
    //     //----------------------------------------------------------------------------------
    // }

    return 0;
}
