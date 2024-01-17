#include <cassert>
#include <cstddef>

#include "fem_unit.hpp"
#include "mandos.hpp"
#include "memory_pool.hpp"
#include "mesh.hpp"
#include "raylib.h"
#include "raylib_imgui.hpp"
#include "rlgl.h"
#include "spring.hpp"
#include "utility_functions.hpp"
#include "viewmandos.hpp"

#define RB_COLOR YELLOW
#define FEM_COLOR (Color){ 255, 109, 194, 255 }
#define PARTICLE_COLOR BLUE
#define MASS_SPRING_COLOR GREEN
#define FROZEN_PARTICLE_COLOR WHITE

Camera3D create_camera(unsigned int FPS = 200) {
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 5.0f, 20.0f };  // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type    return camera;

    return camera;
}

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
static Shader normals_shader;
static Shader diffuse_shader;

Material createMaterialFromShader(Shader shader) {
    Material material = { 0 };
    material.maps = (MaterialMap *)calloc(12, sizeof(MaterialMap));

    // Using rlgl default shader
    material.shader = shader;

    // Using rlgl default texture (1x1 pixel, UNCOMPRESSED_R8G8B8A8, 1 mipmap)
    material.maps[MATERIAL_MAP_DIFFUSE].texture = (Texture2D){ rlGetTextureIdDefault(), 1, 1, 1, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8 };

    material.maps[MATERIAL_MAP_DIFFUSE].color = WHITE;    // Diffuse color
    material.maps[MATERIAL_MAP_SPECULAR].color = WHITE;   // Specular color

    return material;
}

MandosViewer::MandosViewer() {
    InitWindow(initialScreenWidth, initialScreenHeight, "Mandos");
    ClearWindowState(FLAG_WINDOW_RESIZABLE);
    camera = create_camera();
    sphere_model = LoadModelFromMesh(GenMeshSphere(0.1, SPHERE_SUBDIVISIONS, SPHERE_SUBDIVISIONS));
    SetTargetFPS(200);

    // Load basic lighting shader
    base_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/lighting.fs");
    normals_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/normals.fs");
    diffuse_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/picker.fs");
    // Get some required shader locations
    base_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(base_shader, "viewPos");
    base_material = createMaterialFromShader(base_shader);
    sphere_model.materialCount = 1;
    sphere_model.materials[0] = base_material;

}

MandosViewer::~MandosViewer() {
    UnloadModel(sphere_model); // Unloads also the material maps
    UnloadShader(base_shader);
    UnloadShader(normals_shader);
    UnloadShader(diffuse_shader);
    CloseWindow();
}

bool MandosViewer::window_should_close() {
    return WindowShouldClose();
}

void MandosViewer::begin_drawing() {
    mem_pool.reset();
    BeginDrawing();
    ClearBackground(RAYWHITE);
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
    // UpdateCamera(&camera, CAMERA_FREE);

    // Update the shader with the camera view vector (points towards { 0.0f, 0.0f, 0.0f })
    float cameraPos[3] = { camera.position.x, camera.position.y, camera.position.z };
    SetShaderValue(base_shader, base_shader.locs[SHADER_LOC_VECTOR_VIEW], cameraPos, SHADER_UNIFORM_VEC3);
}

bool MandosViewer::is_key_pressed(int Key) { return IsKeyPressed(Key); }

void MandosViewer::draw_particle(const ParticleHandle& particle, const PhysicsState& state) {
    const Vec3 position = particle.particle.get_position(state.x);
    DrawModel(sphere_model, vector3_eigen_to_raylib(position), 1, PARTICLE_COLOR);
}

void draw_mesh_color(const Mat4& transform, const MeshGPU& mesh, MemoryPool& mem_pool, Color color) {
    base_material.maps[MATERIAL_MAP_DIFFUSE].color = color;
    Mesh raymesh = MeshGPUtoRaymesh(mesh, mem_pool);
    DrawMesh(raymesh, base_material, matrix_eigen_to_raylib(transform));
    base_material.maps[MATERIAL_MAP_DIFFUSE].color = WHITE;
}

void MandosViewer::draw_rigid_body(const RigidBodyHandle& rb, const PhysicsState& state, const MeshGPU& mesh) {
    draw_mesh_color(rb.get_transformation_matrix(state), mesh, mem_pool, RB_COLOR);
}
void MandosViewer::draw_mesh(const Mat4& transform, const MeshGPU& mesh) {
    draw_mesh_color(transform, mesh, mem_pool, WHITE);
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
    for (size_t i = 0; i < simulation.simulables.particles.size(); i++) {
        Particle particle = simulation.simulables.particles[i];
        Vec3 position = particle.get_position(state.x);
        Color color = PARTICLE_COLOR;
        for (unsigned int i = 0; i < simulation.frozen_dof.size(); i++) {
            unsigned int frozen_index = simulation.frozen_dof[i];
            if (frozen_index == particle.index) {
                color = WHITE;
                break;
            }
        }
        DrawModel(sphere_model, vector3_eigen_to_raylib(position), 1, color);
    }
}

void MandosViewer::draw_FEM(const FEMHandle& fem, const PhysicsState& state, MeshGPU& gpuMesh, RenderMesh& renderMesh, SimulationMesh& simMesh) {
    const unsigned int dof_index = fem.bounds.dof_index;
    const unsigned int nDoF = fem.bounds.nDoF;

    assert(simMesh.vertices.size() <= nDoF);

    // We should update the simulation mesh from the tetrahedra
    for (unsigned int i = 0; i < simMesh.vertices.size(); i++) {
        simMesh.vertices[i] = state.x[dof_index+i];
    }

    renderMesh.updateFromSimulationMesh(simMesh);
    gpuMesh.updateData(renderMesh);
    draw_mesh_color(Mat4::Identity(), gpuMesh, mem_pool, FEM_COLOR);
}

void MandosViewer::draw_MassSpring(const MassSpringHandle& mass_spring, const PhysicsState& state, MeshGPU& gpuMesh, RenderMesh& renderMesh, SimulationMesh& simMesh) {
    const unsigned int dof_index = mass_spring.bounds.dof_index;
    const unsigned int nDoF = mass_spring.bounds.nDoF;

    assert(simMesh.vertices.size() == nDoF);

    // We should update the simulation mesh from the tetrahedra
    for (unsigned int i = 0; i < simMesh.vertices.size(); i++) {
        simMesh.vertices[i] = state.x[dof_index+i];
    }

    renderMesh.updateFromSimulationMesh(simMesh);
    gpuMesh.updateData(renderMesh);
    draw_mesh_color(Mat4::Identity(), gpuMesh, mem_pool, MASS_SPRING_COLOR);

}


void MandosViewer::draw_FEM_tetrahedrons(const Simulation& simulation, const PhysicsState& state) {
#define MAT(type, name)                                                 \
    for (unsigned int i = 0; i < simulation.energies.fem_elements_##name.size(); i++) { \
        const FEM_Element3D<type>& e = simulation.energies.fem_elements_##name[i]; \
        const Vec3& x1 = e.p1.get_position(state.x);                    \
        const Vec3& x2 = e.p2.get_position(state.x);                    \
        const Vec3& x3 = e.p3.get_position(state.x);                    \
        const Vec3& x4 = e.p4.get_position(state.x);                    \
                                                                        \
        DrawLine3D(vector3_eigen_to_raylib(x1), vector3_eigen_to_raylib(x2), BLUE); \
        DrawLine3D(vector3_eigen_to_raylib(x1), vector3_eigen_to_raylib(x3), BLUE); \
        DrawLine3D(vector3_eigen_to_raylib(x1), vector3_eigen_to_raylib(x4), BLUE); \
        DrawLine3D(vector3_eigen_to_raylib(x4), vector3_eigen_to_raylib(x2), BLUE); \
        DrawLine3D(vector3_eigen_to_raylib(x4), vector3_eigen_to_raylib(x3), BLUE); \
        DrawLine3D(vector3_eigen_to_raylib(x2), vector3_eigen_to_raylib(x3), BLUE); \
    }
    FEM_MATERIAL_MEMBERS
#undef MAT
}
