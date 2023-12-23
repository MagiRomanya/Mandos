#include "memory_pool.hpp"
#include "mesh.hpp"
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include "utility_functions.hpp"
#include "viewmandos.hpp"
#include "render/simulation_visualization.hpp"


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
    std::vector<unsigned short> indices = std::vector<unsigned short>(render_mesh.indices.begin(), render_mesh.indices.end());
    Mesh mesh = {0};
    mesh.vertexCount = render_mesh.vertices.size() / 3;
    mesh.triangleCount = indices.size() / 3;
    mesh.vertices = (float *)std::memcpy(pool.allocate(render_mesh.vertices.size() * sizeof(float)), render_mesh.vertices.data(), render_mesh.vertices.size() * sizeof(float));
    mesh.texcoords = (float *)std::memcpy(pool.allocate(render_mesh.texcoord.size() * sizeof(float)), render_mesh.texcoord.data(), render_mesh.texcoord.size() * sizeof(float));
    mesh.normals = (float *)std::memcpy(pool.allocate(render_mesh.normals.size() * sizeof(float)), render_mesh.normals.data(), render_mesh.normals.size() * sizeof(float));
    mesh.indices = (unsigned short *)std::memcpy(pool.allocate(indices.size() * sizeof(unsigned short)), indices.data(), indices.size() * sizeof(unsigned short));
    UploadMesh(&mesh, false);
    return mesh;
}

void UnloadGPUMesh(const Mesh& mesh) {

    // Unload rlgl mesh vboId data
    rlUnloadVertexArray(mesh.vaoId);

    #define MAX_MESH_VERTEX_BUFFERS 7
    if (mesh.vboId != NULL) for (int i = 0; i < MAX_MESH_VERTEX_BUFFERS; i++) rlUnloadVertexBuffer(mesh.vboId[i]);
    RL_FREE(mesh.vboId);
}

static Camera3D camera;

static Model sphere_model;
static Material base_material;
static Shader base_shader;

MandosViewer::MandosViewer() {
    InitWindow(initialScreenWidth, initialScreenHeight, "Mandos");
    camera = create_camera();
    sphere_model = LoadModelFromMesh(GenMeshSphere(0.1, 15, 15));
    SetTargetFPS(200);

    // Load basic lighting shader
    base_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/lighting.fs");
    // Get some required shader locations
    base_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(base_shader, "viewPos");
    base_material = LoadMaterialDefault();
    base_material.shader = base_shader;
    sphere_model.materialCount = 1;
    sphere_model.materials[0] = base_material;
    SetTraceLogLevel(TraceLogLevel::LOG_WARNING);
}

MandosViewer::~MandosViewer() {
    UnloadModel(sphere_model);
    UnloadMaterial(base_material);
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

void MandosViewer::begin_ImGUI_mode() { begin_ImGUI_mode(); }

void MandosViewer::end_ImGUI_mode() { end_ImGUI_mode(); }

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

void MandosViewer::draw_rigid_body(const RigidBodyHandle& rb, const PhysicsState& state, const RenderMesh& mesh) {
    const Mesh raymesh = RenderMesh_to_RaylibMesh(mesh, mem_pool);
    const Matrix transformation = matrix_eigen_to_raylib(rb.get_transformation_matrix(state));
    Material material = base_material;
    material.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    DrawMesh(raymesh, material, transformation);
    UnloadGPUMesh(raymesh);
}

void MandosViewer::draw_rigid_body(const RigidBodyHandle& rb, const PhysicsState& state, const SimulationMesh& mesh) {
    const Mesh raymesh = SimulationMesh_to_RaylibMesh(mesh, mem_pool);
    const Matrix transformation = matrix_eigen_to_raylib(rb.get_transformation_matrix(state));
    Material material = base_material;
    material.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    DrawMesh(raymesh, material, transformation);
    UnloadGPUMesh(raymesh);
}

void MandosViewer::draw_mesh(const Mat4& transform, const SimulationMesh& mesh) {
    const Mesh raymesh = SimulationMesh_to_RaylibMesh(mesh, mem_pool);
    Material material = LoadMaterialDefault();
    material.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    DrawMesh(raymesh, material, matrix_eigen_to_raylib(transform));
    UnloadGPUMesh(raymesh);
}

void MandosViewer::draw_mesh(const Mat4& transform, const RenderMesh& mesh) {
    const Mesh raymesh = RenderMesh_to_RaylibMesh(mesh, mem_pool);
    Material material = LoadMaterialDefault();
    material.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    DrawMesh(raymesh, material, matrix_eigen_to_raylib(transform));
    UnloadGPUMesh(raymesh);
}

void MandosViewer::draw_mesh(SimulableBounds& bounds, const PhysicsState& state, const SimulationMesh& mesh) {
    static_assert(std::is_same<Scalar, float>::value, "We need scalars to be floats in this implementation.");
    Mesh raymesh = SimulationMesh_to_RaylibMesh(mesh, mem_pool);
    const void* vertices_start = (void*)(state.x.data() + bounds.dof_index);
    memcpy(raymesh.vertices, vertices_start, sizeof(float)*bounds.nDoF);

    Material material = LoadMaterialDefault();
    material.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    DrawMesh(raymesh, material, MatrixIdentity());
    UnloadGPUMesh(raymesh);
}
