#include "mandos.hpp"
#include "viewmandos.hpp"
#include "../mesh.hpp"
#include "../rigid_body.hpp"
#include "../simulation.hpp"
#include "../physics_state.hpp"
#include "../differentiable.hpp"
#include "../async_simulation_loop.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pymandos, m) {
    // BASIC STRUCTURES
    // -----------------------------------------------------------------------------
    py::class_<Simulation>(m, "Simulation")
        .def(py::init())
        .def("get_nDoF", &Simulation::get_nDoF)
        .def("copy",  [](const Simulation &self) {
            return Simulation(self);
        })
        .def_readwrite("initial_state", &Simulation::initial_state)
        .def_readwrite("TimeStep", &Simulation::TimeStep)
        .def_readwrite("MaxNewtonIterations", &Simulation::MaxNewtonIterations)
        .def_readwrite("enable_line_search", &Simulation::enable_line_search)
        ;

    py::class_<PhysicsState>(m, "PhysicsState")
        .def(py::init())
        .def("get_nDoF", &PhysicsState::get_nDoF)
        .def("copy",  [](const PhysicsState &self) {
            return PhysicsState(self);
        })
        .def_readwrite("x", &PhysicsState::x)
        .def_readwrite("v", &PhysicsState::v)
        .def_readwrite("time", &PhysicsState::time)
        ;

    py::class_<EnergyAndDerivatives>(m, "EnergyAndDerivatives")
        .def(py::init<unsigned int>())
        .def_readonly("energy", &EnergyAndDerivatives::energy)
        .def_readonly("gradient", &EnergyAndDerivatives::gradient)
        .def("get_hessian", [](const EnergyAndDerivatives& self) {
            const unsigned int nDoF = self.gradient.size();
            SparseMat hessian(nDoF, nDoF);
            hessian.setFromTriplets(self.hessian_triplets.begin(), self.hessian_triplets.end());
            return hessian;
        })
        ;

    py::class_<RenderMesh>(m, "RenderMesh")
        .def(py::init())
        .def(py::init<std::string, bool>(), py::arg("filename"), py::arg("center") = true)
        .def("updateFromSimulationMesh", &RenderMesh::updateFromSimulationMesh)
        .def("smoothNormals", &RenderMesh::smoothNormals)
        .def_readonly("vertices", &RenderMesh::vertices)
        .def_readonly("texcoords", &RenderMesh::texcoords)
        .def_readonly("normals", &RenderMesh::normals)
        .def_readonly("tangents", &RenderMesh::tangents)
        ;

    py::class_<SimulationMesh>(m, "SimulationMesh")
        .def(py::init())
        .def(py::init<std::string, bool>(), py::arg("filename"), py::arg("center") = true)
        .def(py::init<const RenderMesh&, bool>(), py::arg("simMesh"), py::arg("center") = false)
        .def_readonly("vertices", &SimulationMesh::vertices)
        .def_readonly("indices", &SimulationMesh::indices)
        ;

    py::class_<TetrahedronMesh>(m, "TetrahedronMesh")
        .def(py::init())
        .def_readonly("vertices", &TetrahedronMesh::vertices)
        .def_readonly("indices", &TetrahedronMesh::indices)
        ;

    m.def("load_curve", &LoadCurveTinyObj);

    // SIMULATION
    // -----------------------------------------------------------------------------
    m.def("simulation_step", py::overload_cast<const Simulation&, PhysicsState&, EnergyAndDerivatives&>(&simulation_step));

    m.def("compute_inertia_tensor_PARTICLES", &compute_initial_inertia_tensor_PARTICLES);

    m.def("compute_COM_position_PARTICLES", &compute_COM_position_PARTICLES);

    m.def("join_rigid_body_with_particle", &join_rigid_body_with_particle);

    m.def("join_particles_with_spring", &join_particles_with_spring);

    m.def("join_rigid_body_com_with_spring", &join_rigid_body_com_with_spring);

    m.def("join_rigid_body_with_spring", py::overload_cast<Simulation&, const RigidBodyHandle&, const Vec3&, const RigidBodyHandle&, const Vec3&, Scalar, Scalar, Scalar>(&join_rigid_body_with_spring));

    m.def("join_rigid_body_with_spring", py::overload_cast<Simulation&, const RigidBodyHandle&, const Vec3&, const RigidBodyHandle&, const Vec3&, Scalar, Scalar>(&join_rigid_body_with_spring));

    m.def("join_rigid_body_with_rod_segment", &join_rigid_body_with_rod_segment);

    m.def("get_particle", &get_particle_handle);

    m.def("compute_tetrahedrons", py::overload_cast<const std::vector<unsigned int>&, const std::vector<Scalar>&, TetrahedronMesh&>(&tetgen_compute_tetrahedrons));

    m.def("compute_energy_and_derivatives", [](const Simulation& simulation, const PhysicsState& state, const PhysicsState& state0){
        EnergyAndDerivatives f(state.get_nDoF());
        compute_energy_and_derivatives(simulation.TimeStep, simulation.energies, state, state0, f);
        return f;
    });

    m.def("compute_energy_and_derivatives_finite", [](const Simulation& simulation, const PhysicsState& state, const PhysicsState& state0){
        EnergyAndDerivatives f(state.get_nDoF());
        compute_energy_and_derivatives_finite(simulation.TimeStep, simulation.energies, state, state0, f);
        return f;
    });

    // ASYNC SIMULATION
    // -------------------------------------------------------------------------------------
    m.def("simulation_async_loop", &simulation_async_loop);

    m.def("simulation_async_loop_request_iteration", &simulation_async_loop_request_iteration);

    m.def("get_current_physics_state", &get_current_physics_state);

    m.def("set_current_physics_state", &set_current_physics_state);

    m.def("get_current_energy_and_derivatives", &get_current_energy_and_derivatives);

    // SIMULABLES
    // ----------------------------------------------------------------

    py::class_<RigidBodyHandle>(m, "RigidBody")
        .def(py::init<Simulation&, Scalar, std::vector<Scalar>, bool>(),
             py::arg("simulation"), py::arg("mass"), py::arg("vertices"), py::arg("global") = false)
        .def(py::init<Simulation&, Scalar, Mat3, bool>(),
             py::arg("simulation"), py::arg("mass"), py::arg("inertia_tensor"), py::arg("global") = false)
        .def("set_com_position", &RigidBodyHandle::set_COM_position)
        .def("get_com_position", &RigidBodyHandle::get_COM_position)
        .def("set_com_initial_position", &RigidBodyHandle::set_COM_initial_position)
        .def("set_com_initial_velocity", &RigidBodyHandle::set_COM_initial_velocity)
        .def("set_initial_orientation", &RigidBodyHandle::set_initial_orientation)
        .def("set_initial_angular_velocity", &RigidBodyHandle::set_initial_angular_velocity)
        .def("add_gravity", &RigidBodyHandle::add_gravity)
        .def("freeze_translation", &RigidBodyHandle::freeze_translation)
        .def("freeze_rotation", &RigidBodyHandle::freeze_rotation)
        .def("get_transformation_matrix", &RigidBodyHandle::get_transformation_matrix)
        .def("compute_energy", &RigidBodyHandle::compute_energy)
        ;

    py::class_<MassSpringHandle>(m, "MassSpring")
        .def(py::init<Simulation&, const std::vector<Scalar>&, const std::vector<unsigned int>&, Scalar, Scalar, Scalar, Scalar>())
        .def("compute_center_of_mass", &MassSpringHandle::compute_center_of_mass)
        .def("set_initial_com_position", &MassSpringHandle::set_initial_COM_position)
        .def("freeze_particles", &MassSpringHandle::freeze_particles)
        .def("add_gravity", &MassSpringHandle::add_gravity)
        .def("get_n_particles", &MassSpringHandle::get_n_particles)
        .def("get_dof_vector", &MassSpringHandle::get_dof_vector)
        ;

    py::class_<ParticleHandle>(m, "Particle")
        .def(py::init<Simulation&, Scalar>())
        .def("set_initial_position", &ParticleHandle::set_initial_position)
        .def("set_initial_velocity", &ParticleHandle::set_initial_velocity)
        .def("add_gravity", &ParticleHandle::add_gravity)
        .def("freeze", &ParticleHandle::freeze)
        .def("get_position", &ParticleHandle::get_position)
        ;

    py::class_<FEMHandle>(m, "FEM")
        .def(py::init<Simulation&, const std::vector<Scalar>, const std::vector<unsigned int>, Scalar, Scalar, Scalar>())
        .def("compute_center_of_mass", &FEMHandle::compute_center_of_mass)
        .def("freeze_particles", &FEMHandle::freeze_particles)
        .def("get_n_particles", &FEMHandle::get_n_particles)
        .def("add_gravity", &FEMHandle::add_gravity)
        .def("set_com_position", &FEMHandle::set_com_position)
        ;

    py::class_<RodSegmentParameters>(m, "RodParameters")
        .def(py::init())
        .def_readwrite("Ks", &RodSegmentParameters::Ks)
        .def_readwrite("translational_damping", &RodSegmentParameters::translational_damping)
        .def_readwrite("constraint_stiffness", &RodSegmentParameters::constraint_stiffness)
        .def_readwrite("intrinsic_darboux", &RodSegmentParameters::intrinsic_darboux)
        .def_readwrite("stiffness_tensor", &RodSegmentParameters::stiffness_tensor)
        ;

    py::class_<RodHandle>(m, "Rod")
        .def(py::init<Simulation&, unsigned int, const Vec3&, Scalar, RodSegmentParameters>())
        .def(py::init<Simulation&, const std::vector<Scalar>&, Scalar, RodSegmentParameters>())
        .def("compute_center_of_mass", &RodHandle::compute_center_of_mass)
        .def("add_gravity", &RodHandle::add_gravity)
        .def("set_rigid_body_initial_position", &RodHandle::set_rigid_body_initial_position)
        .def("set_rigid_body_initial_velocity", &RodHandle::set_rigid_body_initial_velocity)
        .def("set_rigid_body_initial_orientation", &RodHandle::set_rigid_body_initial_orientation)
        .def("set_rigid_body_initial_angular_velocity", &RodHandle::set_rigid_body_initial_angular_velocity)
        .def("set_initial_rod_position", &RodHandle::set_initial_rod_position)
        .def("get_rigid_body", &RodHandle::get_rigid_body)
        .def("freeze_rigid_body", &RodHandle::freeze_rigid_body)
        ;

    py::class_<RigidBodySpringHandle>(m,"RigidBodySpring")
        .def(py::init<Simulation&, const RigidBodyHandle&, const Vec3&, const RigidBodyHandle&, const Vec3&, Scalar, Scalar, Scalar>())
        .def(py::init<Simulation&, const RigidBodyHandle&, const Vec3&, const RigidBodyHandle&, const Vec3&, Scalar, Scalar>())
        .def("set_rest_length", &RigidBodySpringHandle::set_rest_length)
        .def("set_stiffness", &RigidBodySpringHandle::set_stiffness)
        ;

    // COLLIDERS
    // -----------------------------------------------------------------------------
    py::class_<PlaneColliderHandle>(m, "PlaneCollider")
        .def(py::init<Simulation&>())
        .def("set_origin_position", &PlaneColliderHandle::set_origin_position)
        .def("set_direction", &PlaneColliderHandle::set_direction)
        ;

    py::class_<SphereColliderHandle>(m, "SphereCollider")
        .def(py::init<Simulation&>())
        .def("set_origin_position", &SphereColliderHandle::set_origin_position)
        .def("set_radius", &SphereColliderHandle::set_radius)
        ;

    py::class_<SDFColliderHandle>(m, "SDFCollider")
        .def(py::init<Simulation&, const SimulationMesh&>())
        ;

    // DIFFERENTIABLE SIMULATION
    // -----------------------------------------------------------------------------
    py::class_<LossFunctionAndDerivatives>(m, "LossFunctionAndDerivatives")
        .def(py::init<>())
        .def_readwrite("loss", &LossFunctionAndDerivatives::loss)
        .def_readwrite("loss_position_partial_derivative", &LossFunctionAndDerivatives::loss_position_partial_derivative)
        .def_readwrite("loss_velocity_partial_derivative", &LossFunctionAndDerivatives::loss_velocity_partial_derivative)
        .def_readwrite("loss_parameter_partial_derivative", &LossFunctionAndDerivatives::loss_parameter_partial_derivative)
        ;

    m.def("compute_loss_function_gradient_backpropagation", &compute_loss_function_gradient_backpropagation,
          "A function that computes the loss function gradient with respect to parameters using back propagation of gradients.",
          py::arg("simulation"), py::arg("trajectory"), py::arg("loss"), py::arg("dx0_dp"), py::arg("dv0_dp"), py::arg("") = 0);

    m.def("compute_loss_function_gradient_backpropagation_control", &compute_loss_function_gradient_backpropagation_control,
          "A function that computes the loss function gradient with respect to parameters using back propagation of gradients.",
          py::arg("simulation"), py::arg("trajectory"), py::arg("loss"), py::arg("") = 0);

    m.def("compute_loss_function_gradient_backpropagation_1_step_position", &compute_loss_function_gradient_backpropagation_1_step_position);

    m.def("compute_loss_function_gradient_backpropagation_1_step_velocity", &compute_loss_function_gradient_backpropagation_1_step_velocity);


    // RENDERING
    // -----------------------------------------------------------------------------
#ifdef USE_VISUALIZER
    py::class_<MeshGPU>(m, "MeshGPU")
        .def(py::init<const RenderMesh&>())
        .def("updateData", &MeshGPU::updateData)
        ;
    py::class_<SkinnedRodGPU>(m, "SkinnedRodGPU")
        .def(py::init<const RenderMesh&, const unsigned int>())
        .def(py::init<const RenderMesh&, const RodHandle&>())
        ;

    py::class_<MandosViewer>(m, "MandosViewer")
        .def(py::init<>())
        .def("disable_render_logs", &MandosViewer::disable_render_logs)
        .def("window_should_close", &MandosViewer::window_should_close)
        .def("begin_drawing", &MandosViewer::begin_drawing)
        .def("end_drawing", &MandosViewer::end_drawing)
        .def("is_key_pressed", &MandosViewer::is_key_pressed)
        .def("draw_particle", &MandosViewer::draw_particle)
        .def("draw_rigid_body", &MandosViewer::draw_rigid_body)
        .def("draw_FEM", &MandosViewer::draw_FEM)
        .def("draw_MassSpring", &MandosViewer::draw_MassSpring)
        .def("draw_rod", &MandosViewer::draw_rod)
        .def("draw_mesh", &MandosViewer::draw_mesh)
        .def("draw_springs", &MandosViewer::draw_springs)
        .def("draw_rods", &MandosViewer::draw_rods)
        .def("draw_particles", &MandosViewer::draw_particles)
        .def("draw_rigid_bodies", &MandosViewer::draw_rigid_bodies)
        .def("draw_FEM_tetrahedrons", &MandosViewer::draw_FEM_tetrahedrons_lines)
        // .def("draw_particle_indices", &MandosViewer::draw_particle_indices)
        .def("draw_simulation_state", &MandosViewer::draw_simulation_state)
        .def("draw_vector", &MandosViewer::draw_vector)
        ;

#endif // USE_VISUALIZER
}
