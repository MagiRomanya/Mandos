set(HEADERS
  async_simulation_loop.hpp
  colliders.hpp
  differentiable.hpp
  edge.hpp
  energies.hpp
  fem_element.hpp
  gravity.hpp
  inertia_energies.hpp
  integrators.hpp
  linear_algebra.hpp
  mesh.hpp
  particle.hpp
  particle_rigid_body_coupling.hpp
  physics_state.hpp
  rigid_body.hpp
  rod_segment.hpp
  simulable_generator.hpp
  simulation.hpp
  sinc.hpp
  spring.hpp
  utility_functions.hpp
  mandos.hpp
)

set(SOURCES
  async_simulation_loop.cpp
  colliders.cpp
#   compute_tetrahedrons.cpp  
  differentiable.cpp
  edge.cpp
  fem_element.cpp
  inertia_energies.cpp
  integrators.cpp
  linear_algebra.cpp
  mandos.cpp
  particle_rigidbody_coupling.cpp
  rigid_body.cpp
  rod_segment.cpp
  simulation.cpp
  sinc.cpp
  spring.cpp
)

add_library(Mandos SHARED ${HEADERS} ${SOURCES})

# Add the src directory as the include directory. Includes should be added as <Mandos/Foo.hpp>
target_include_directories(Mandos PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>
)

target_link_libraries(Mandos 
PRIVATE 
  Eigen3::Eigen
  spdlog::spdlog
)

generate_export_header(Mandos)
