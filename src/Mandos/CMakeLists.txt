set(HEADERS
  async_simulation_loop.hpp
#   colliders.hpp
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
  particle_rigid_body_copuling.hpp
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
#   colliders.cpp
#   compute_tetrahedrons.cpp  
  differentiable.cpp
  edge.cpp
  fem_element.cpp
  inertia_energies.cpp
  integrators.cpp
  mandos.cpp
  mesh_boundary.cpp
  particle_rigidbody_copuling.cpp
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
)

generate_export_header(Mandos)
