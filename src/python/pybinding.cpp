#include "mandos.hpp"
#include "mesh.hpp"
#include "rigid_body.hpp"
#include "simulation.hpp"
#include "viewmandos.hpp"
#include "physics_state.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pymandos, m) {
    // BASIC STRUCTURES
    // -----------------------------------------------------------------------------
    py::class_<Simulation>(m, "Simulation")
        .def(py::init())
        .def_readonly("initial_state", &Simulation::initial_state)
        ;

    py::class_<PhysicsState>(m, "PhysicsState")
        .def(py::init())
        .def_readonly("x", &PhysicsState::x)
        .def_readonly("x_old", &PhysicsState::x_old)
        .def("copy",  [](const PhysicsState &self) {
            return PhysicsState(self);
        })
        ;

    py::class_<EnergyAndDerivatives>(m, "EnergyAndDerivatives")
        .def(py::init<unsigned int>())
        .def_readonly("energy", &EnergyAndDerivatives::energy)
        .def_readonly("gradient", &EnergyAndDerivatives::gradient)
        ;

    py::class_<RenderMesh>(m, "RenderMesh")
        .def(py::init())
        .def(py::init<std::string>())
        .def("updateFromSimulationMesh", &RenderMesh::updateFromSimulationMesh)
        .def_readonly("vertices", &RenderMesh::vertices)
        .def_readonly("texcoord", &RenderMesh::texcoord)
        .def_readonly("normals", &RenderMesh::normals)
        ;

    py::class_<SimulationMesh>(m, "SimulationMesh")
        .def(py::init())
        .def(py::init<std::string>())
        .def(py::init<const RenderMesh&>())
        .def_readonly("vertices", &SimulationMesh::vertices)
        .def_readonly("indices", &SimulationMesh::indices)
        ;

    py::class_<TetrahedronMesh>(m, "TetrahedronMesh")
        .def(py::init())
        .def_readonly("vertices", &TetrahedronMesh::vertices)
        .def_readonly("indices", &TetrahedronMesh::indices)
        ;


    // SIMULATION
    // -----------------------------------------------------------------------------
    m.def("simulation_step", py::overload_cast<const Simulation&, PhysicsState&, EnergyAndDerivatives&>(&simulation_step));

    m.def("compute_inertia_tensor_PARTICLES", &compute_initial_inertia_tensor_PARTICLES);

    m.def("compute_COM_position_PARTICLES", &compute_COM_position_PARTICLES);

    m.def("join_rigid_body_with_particle", &join_rigid_body_with_particle);

    m.def("join_particles_with_spring", &join_particles_with_spring);

    m.def("compute_tetrahedrons", py::overload_cast<const std::vector<unsigned int>&, const std::vector<float>&, TetrahedronMesh&>(&tetgen_compute_tetrahedrons));

    py::class_<RigidBodyHandle>(m, "RigidBody")
        .def(py::init<Simulation&, Scalar, std::vector<Scalar>>())
        .def(py::init<Simulation&, Scalar, Mat3>())
        .def("set_com_initial_position", &RigidBodyHandle::set_COM_initial_position)
        .def("set_com_initial_velocity", &RigidBodyHandle::set_COM_initial_velocity)
        .def("set_initial_orientation", &RigidBodyHandle::set_initial_orientation)
        .def("set_initial_angular_velocity", &RigidBodyHandle::set_initial_angular_velocity)
        .def("add_gravity", &RigidBodyHandle::add_gravity)
        .def("freeze_translation", &RigidBodyHandle::freeze_translation)
        .def("freeze_rotation", &RigidBodyHandle::freeze_rotation)
        .def("get_transformation_matrix", &RigidBodyHandle::get_transformation_matrix)
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
        ;

    // RENDERING
    // -----------------------------------------------------------------------------
    py::class_<MeshGPU>(m, "MeshGPU")
        .def(py::init<const RenderMesh&>())
        .def("updateData", &MeshGPU::updateData)
        ;

    py::class_<MandosViewer>(m, "MandosViewer")
        .def(py::init<>())
        .def("window_should_close", &MandosViewer::window_should_close)
        .def("begin_drawing", &MandosViewer::begin_drawing)
        .def("end_drawing", &MandosViewer::end_drawing)
        .def("begin_3D_mode", &MandosViewer::begin_3D_mode)
        .def("end_3D_mode", &MandosViewer::end_3D_mode)
        .def("update_camera", &MandosViewer::update_camera)
        .def("is_key_pressed", &MandosViewer::is_key_pressed)
        .def("draw_particle", &MandosViewer::draw_particle)
        .def("draw_rigid_body", &MandosViewer::draw_rigid_body)
        .def("draw_FEM", &MandosViewer::draw_FEM)
        .def("draw_MassSpring", &MandosViewer::draw_MassSpring)
        .def("draw_mesh", &MandosViewer::draw_mesh)
        .def("draw_springs", &MandosViewer::draw_springs)
        .def("draw_particles", &MandosViewer::draw_particles)
        .def("draw_FEM_tetrahedrons", &MandosViewer::draw_FEM_tetrahedrons)
        ;

}
