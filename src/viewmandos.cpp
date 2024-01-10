#include "mandos.hpp"
#include "memory_pool.hpp"
#include "mesh.hpp"
#include "raylib.h"
#include "raylib_imgui.hpp"
#include "rlgl.h"
#include "spring.hpp"
#include "utility_functions.hpp"
#include "viewmandos.hpp"
#include "render/simulation_visualization.hpp"
#include <cassert>
#include <cstddef>
#include <cstdlib>

Mesh SimulationMesh_to_RaylibMesh(const SimulationMesh& sim_mesh, MemoryPool& pool) {
    std::vector<unsigned short> indices = std::vector<unsigned short>(sim_mesh.indices.begin(), sim_mesh.indices.end());
    Mesh mesh = {0};
    mesh.vertexCount = sim_mesh.vertices.size() / 3;
    mesh.triangleCount = indices.size() / 3;
    mesh.vertices = (float *)std::memcpy(pool.allocate(sim_mesh.vertices.size() * sizeof(float)), sim_mesh.vertices.data(), sim_mesh.vertices.size() * sizeof(float));
    mesh.texcoords = NULL;
    mesh.normals = NULL;
    mesh.indices = (unsigned short *)std::memcpy(pool.allocate(indices.size() * sizeof(unsigned short)), indices.data(), indices.size() * sizeof(unsigned short));
    UploadMesh(&mesh, false);
    return mesh;
}

Mesh RenderMesh_to_RaylibMesh(const RenderMesh& render_mesh, MemoryPool& pool) {
    Mesh mesh = {0};
    mesh.vertexCount = render_mesh.vertices.size() / 3; // 3 coordinates per vertex
    mesh.triangleCount = mesh.vertexCount / 3; // 3 vertices per triangle
    mesh.vertices = (float *)std::memcpy(pool.allocate(render_mesh.vertices.size() * sizeof(float)), render_mesh.vertices.data(), render_mesh.vertices.size() * sizeof(float));
    mesh.texcoords = (float *)std::memcpy(pool.allocate(render_mesh.texcoord.size() * sizeof(float)), render_mesh.texcoord.data(), render_mesh.texcoord.size() * sizeof(float));
    mesh.normals = (float *)std::memcpy(pool.allocate(render_mesh.normals.size() * sizeof(float)), render_mesh.normals.data(), render_mesh.normals.size() * sizeof(float));
    mesh.indices = NULL;
    UploadMesh(&mesh, false);
    return mesh;
}

MeshGPU::MeshGPU(const RenderMesh& mesh) {
    // Allocate resources
    vertices = (float*) calloc(sizeof(float), mesh.vertices.size());
    texcoords = (float*) calloc(sizeof(float), mesh.texcoord.size());
    normals = (float*) calloc(sizeof(float), mesh.normals.size());

    // Copy the data
    nVertices = mesh.vertices.size() / 3;
    vertices = (float *) std::memcpy(vertices, mesh.vertices.data(), mesh.vertices.size()*sizeof(float));
    texcoords = (float *) std::memcpy(texcoords, mesh.texcoord.data(), mesh.texcoord.size()*sizeof(float));
    normals = (float *) std::memcpy(normals, mesh.normals.data(), mesh.normals.size()*sizeof(float));
    if (!vertices) std::cerr << "MeshGPU::MeshGPU: vertices null" << std::endl;
    if (!texcoords) std::cerr << "MeshGPU::MeshGPU: texcoord null" << std::endl;
    if (!normals) std::cerr << "MeshGPU::MeshGPU: normals null" << std::endl;

    // Uload the mesh to GPU
    Mesh raymesh = {0};
    raymesh.vertexCount = mesh.vertices.size() / 3; // 3 coordinates per vertex
    raymesh.triangleCount = raymesh.vertexCount / 3; // 3 vertices per triangle
    raymesh.vertices = vertices;
    raymesh.texcoords = texcoords;
    raymesh.normals = normals;
    UploadMesh(&raymesh, true);
    verticesVBO = raymesh.vboId[0];
    texcoordsVBO = raymesh.vboId[1];
    normalsVBO = raymesh.vboId[2];
    VAO = raymesh.vaoId;
    RL_FREE(raymesh.vboId);
}

void MeshGPU::updateData(const RenderMesh& mesh) {
    vertices = (float*) std::memcpy(vertices, mesh.vertices.data(), mesh.vertices.size()*sizeof(float));
    normals = (float*) std::memcpy(normals, mesh.normals.data(), mesh.normals.size()*sizeof(float));
    rlUpdateVertexBuffer(verticesVBO, vertices, mesh.vertices.size() * sizeof(float), 0);
    rlUpdateVertexBuffer(normalsVBO, normals, mesh.normals.size() * sizeof(float), 0);
}

MeshGPU::~MeshGPU() {
    free(vertices);
    free(texcoords);
    free(normals);

    // Unload the mesh buffers
    rlUnloadVertexArray(VAO);
    rlUnloadVertexBuffer(verticesVBO);
    rlUnloadVertexBuffer(texcoordsVBO);
    rlUnloadVertexBuffer(normalsVBO);
}

Mesh MeshGPUtoRaymesh(const MeshGPU& mesh, MemoryPool& pool) {
    Mesh raymesh = {0};
    raymesh.vertexCount = mesh.nVertices;
    raymesh.triangleCount = raymesh.vertexCount / 3; // 3 vertices per triangle
    raymesh.vertices = mesh.vertices;
    raymesh.texcoords = mesh.texcoords;
    raymesh.normals = mesh.normals;
    // Copy the vaos
    raymesh.vaoId = mesh.VAO;
    raymesh.vboId = (unsigned int*) pool.allocate(sizeof(unsigned int) * 3);
    raymesh.vboId[0] = mesh.verticesVBO;
    raymesh.vboId[1] = mesh.texcoordsVBO;
    raymesh.vboId[2] = mesh.normalsVBO;

    return raymesh;
}

void UnloadGPUMesh(const Mesh& mesh) {

    // Unload rlgl mesh vboId data
    rlUnloadVertexArray(mesh.vaoId);

    #define MAX_MESH_VERTEX_BUFFERS 7
    if (mesh.vboId != NULL) for (int i = 0; i < MAX_MESH_VERTEX_BUFFERS; i++) rlUnloadVertexBuffer(mesh.vboId[i]);
    RL_FREE(mesh.vboId);
    #undef MAX_MESH_VERTEX_BUFFERS
}

static Camera3D camera;

#define SPHERE_SUBDIVISIONS 30
static Model sphere_model;
static Material base_material;
static Shader base_shader;

MandosViewer::MandosViewer() {
    InitWindow(initialScreenWidth, initialScreenHeight, "Mandos");
    camera = create_camera();
    sphere_model = LoadModelFromMesh(GenMeshSphere(0.1, SPHERE_SUBDIVISIONS, SPHERE_SUBDIVISIONS));
    SetTargetFPS(200);

    // Load basic lighting shader
    base_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/lighting.fs");
    // Get some required shader locations
    base_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(base_shader, "viewPos");
    base_material = LoadMaterialDefault();
    base_material.shader = base_shader;
    sphere_model.materialCount = 1;
    sphere_model.materials[0] = base_material;
    // SetTraceLogLevel(TraceLogLevel::LOG_WARNING);
}

MandosViewer::~MandosViewer() {
    UnloadModel(sphere_model);
    UnloadShader(base_shader);
    CloseWindow();
}

bool MandosViewer::window_should_close() {
    return WindowShouldClose();
}

void MandosViewer::begin_drawing() {
    mem_pool.reset();
    BeginDrawing();
    ClearBackground(WHITE);
}

void MandosViewer::end_drawing() {
    DrawFPS(10, 10);
    EndDrawing();
}

void MandosViewer::begin_3D_mode() {
    BeginMode3D(camera);
    DrawGrid(30, 1.0f);
}

void MandosViewer::end_3D_mode() { EndMode3D(); }


void MandosViewer::begin_ImGUI_mode() { ImGuiBeginDrawing(); }

void MandosViewer::end_ImGUI_mode() { ImGuiEndDrawing(); }

void MandosViewer::update_camera() {
    UpdateCamera(&camera, CAMERA_FREE);

    // Update the shader with the camera view vector (points towards { 0.0f, 0.0f, 0.0f })
    float cameraPos[3] = { camera.position.x, camera.position.y, camera.position.z };
    SetShaderValue(base_shader, base_shader.locs[SHADER_LOC_VECTOR_VIEW], cameraPos, SHADER_UNIFORM_VEC3);
}

bool MandosViewer::is_key_pressed(int Key) { return IsKeyPressed(Key); }

void MandosViewer::draw_particle(const ParticleHandle& particle, const PhysicsState& state) {
    const Vec3 position = particle.particle.get_position(state.x);
    DrawModel(sphere_model, vector3_eigen_to_raylib(position), 1, DARKPURPLE);
}

void MandosViewer::draw_rigid_body(const RigidBodyHandle& rb, const PhysicsState& state, const MeshGPU& mesh) {
    this->draw_mesh(rb.get_transformation_matrix(state), mesh);
}

void MandosViewer::draw_mesh(const Mat4& transform, const MeshGPU& mesh) {
    Material material = base_material;
    material.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    Mesh raymesh = MeshGPUtoRaymesh(mesh, mem_pool);
    DrawMesh(raymesh, material, matrix_eigen_to_raylib(transform));
}

void MandosViewer::draw_springs(const Simulation& simulation, const PhysicsState& state) {
    for (size_t i = 0; i < simulation.energies.particle_springs.size(); i++) {
        ParticleSpring s = simulation.energies.particle_springs[i];
        Vec3 x1 = s.p1.get_position(state.x);
        Vec3 x2 = s.p2.get_position(state.x);
        DrawLine3D(vector3_eigen_to_raylib(x1), vector3_eigen_to_raylib(x2), BLACK);
    }
}

void MandosViewer::draw_particles(const Simulation& simulation, const PhysicsState& state) {
    DEBUG_LOG(simulation.simulables.particles.size());
    for (size_t i = 0; i < simulation.simulables.particles.size(); i++) {
        Particle particle = simulation.simulables.particles[i];
        Vec3 position = particle.get_position(state.x);
        Color color = GREEN;
        for (unsigned int i = 0; i < simulation.frozen_dof.size(); i++) {
            unsigned int frozen_index = simulation.frozen_dof[i];
            if (frozen_index == particle.index) {
                color = RED;
                break;
            }
        }
        DrawModel(sphere_model, vector3_eigen_to_raylib(position), 1, color);
    }
}

void MandosViewer::draw_FEM(const FEMHandle& fem, const PhysicsState& state, MeshGPU& gpuMesh, RenderMesh& renderMesh, SimulationMesh& simMesh) {
    std::cerr << "ERROR: MandosViewer::draw_FEM: TODO!" << std::endl;
    exit(1);
    const unsigned int dof_index = fem.bounds.dof_index;
    const unsigned int nDoF = fem.bounds.nDoF;

    // We should update the simulation mesh from the tetrahedra

    renderMesh.updateFromSimulationMesh(simMesh);
    gpuMesh.updateData(renderMesh);
}


void MandosViewer::draw_FEM_tetrahedrons(const Simulation& simulation, const PhysicsState& state) {
    for (unsigned int i = 0; i < simulation.energies.fem_elements_3d.size(); i++) {
        const FEM_Element3D& e = simulation.energies.fem_elements_3d[i];
        const Vec3& x1 = e.p1.get_position(state.x);
        const Vec3& x2 = e.p2.get_position(state.x);
        const Vec3& x3 = e.p3.get_position(state.x);
        const Vec3& x4 = e.p4.get_position(state.x);

        DrawLine3D(vector3_eigen_to_raylib(x1), vector3_eigen_to_raylib(x2), BLUE);
        DrawLine3D(vector3_eigen_to_raylib(x1), vector3_eigen_to_raylib(x3), BLUE);
        DrawLine3D(vector3_eigen_to_raylib(x1), vector3_eigen_to_raylib(x4), BLUE);
        DrawLine3D(vector3_eigen_to_raylib(x4), vector3_eigen_to_raylib(x2), BLUE);
        DrawLine3D(vector3_eigen_to_raylib(x4), vector3_eigen_to_raylib(x3), BLUE);
        DrawLine3D(vector3_eigen_to_raylib(x2), vector3_eigen_to_raylib(x3), BLUE);
    }

}
