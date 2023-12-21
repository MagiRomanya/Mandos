#include "mandos.hpp"
#include "physics_state.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_MODULE(pymandos, m) {
    // BASIC STRUCTURES
    // -----------------------------------------------------------------------------
    py::class_<Simulation>(m, "Simulation")
        .def(py::init())
        ;

    py::class_<PhysicsState>(m, "PhysicsState")
        .def(py::init())
        ;

    // EXPOSE SIMULABLES
    // -----------------------------------------------------------------------------
    py::class_<RigidBodyHandle>(m, "RigidBody")
        .def(py::init<Simulation&, Scalar, std::vector<Scalar>>())
        .def(py::init<Simulation&, Scalar, Mat3>())
        .def("set_com_initial_position", &RigidBodyHandle::set_COM_initial_position)
        .def("set_com_initial_velocity", &RigidBodyHandle::set_COM_initial_velocity)
        .def("set_initial_orientation", &RigidBodyHandle::set_initial_orientation)
        .def("set_initial_anglular_velocity", &RigidBodyHandle::set_initial_angluar_velocity)
        .def("add_gravity", &RigidBodyHandle::add_gravity)
        .def("freeze_translation", &RigidBodyHandle::freeze_translation)
        .def("freeze_rotation", &RigidBodyHandle::freeze_rotation)
        .def("get_transformation_matrix", &RigidBodyHandle::get_transformation_matrix)
        ;

    py::class_<MassSpringHandle>(m, "MassSpring")
        .def(py::init<Simulation&, const std::vector<Scalar>&, const std::vector<unsigned int>&, Scalar, Scalar, Scalar, Scalar>())
        .def("compute_center_of_mass", &MassSpringHandle::compute_center_of_mass)
        .def("set_initial_com_position", &MassSpringHandle::set_initial_COM_position)
        .def("freezze_particles", &MassSpringHandle::freeze_particles)
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

    m.def("join_particles_with_spring", &join_particles_with_spring);

    py::class_<FEMHandle>(m, "FEM")
        .def(py::init<Simulation&, const std::vector<Scalar>, const std::vector<unsigned int>, Scalar, Scalar, Scalar>())
        .def("compute_center_of_mass", &FEMHandle::compute_center_of_mass)
        .def("freeze_particles", &FEMHandle::freeze_particles)
        .def("get_n_particles", &FEMHandle::get_n_particles)
        ;
}
