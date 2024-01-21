#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <stdlib.h>
#include <Eigen/Core>
#include <vector>

#include "fem_unit.hpp"
#include "linear_algebra.hpp"
#include "mandos.hpp"
#include "memory_pool.hpp"
#include "mesh.hpp"
#include "raylib.h"
#include "raylib_imgui.hpp"
#include "raymath.h"
#include "rlgl.h"
#include "spring.hpp"
#include "viewmandos.hpp"
#include "utility_functions.hpp"

#define RB_COLOR YELLOW
#define FEM_COLOR PINK
#define PARTICLE_COLOR BLUE
#define MASS_SPRING_COLOR GREEN
#define FROZEN_PARTICLE_COLOR WHITE

inline Matrix matrix_eigen_to_raylib(const Mat4& m) {
    Matrix r = {
    m(0,0), m(0,1), m(0, 2), m(0, 3),
    m(1,0), m(1,1), m(1, 2), m(1, 3),
    m(2,0), m(2,1), m(2, 2), m(2, 3),
    m(3,0), m(3,1), m(3, 2), m(3, 3),
    };
    return r;
}

inline Vector3 vector3_eigen_to_raylib(const Vec3& v) {
    return Vector3{v.x(), v.y(), v.z()};
}

Camera3D create_camera() {
    Camera3D camera = { 0 };
    camera.position = Vector3( 0.0f, 5.0f, 20.0f );  // Camera position
    camera.target = Vector3( 0.0f, 0.0f, 0.0f );      // Camera looking at point
    camera.up = Vector3( 0.0f, 1.0f, 0.0f );          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type    return camera;
    return camera;
}

// DEPRECATED
// Mesh SimulationMesh_to_RaylibMesh(const SimulationMesh& sim_mesh, MemoryPool& pool) {
//     std::vector<unsigned short> indices = std::vector<unsigned short>(sim_mesh.indices.begin(), sim_mesh.indices.end());
//     Mesh mesh = {0};
//     mesh.vertexCount = sim_mesh.vertices.size() / 3;
//     mesh.triangleCount = indices.size() / 3;
//     mesh.vertices = (float *)std::memcpy(pool.allocate(sim_mesh.vertices.size() * sizeof(float)), sim_mesh.vertices.data(), sim_mesh.vertices.size() * sizeof(float));
//     mesh.texcoords = NULL;
//     mesh.normals = NULL;
//     mesh.indices = (unsigned short *)std::memcpy(pool.allocate(indices.size() * sizeof(unsigned short)), indices.data(), indices.size() * sizeof(unsigned short));
//     UploadMesh(&mesh, false);
//     return mesh;
// }

// DEPRECATED
// Mesh RenderMesh_to_RaylibMesh(const RenderMesh& render_mesh, MemoryPool& pool) {
//     Mesh mesh = {0};
//     mesh.vertexCount = render_mesh.vertices.size() / 3; // 3 coordinates per vertex
//     mesh.triangleCount = mesh.vertexCount / 3; // 3 vertices per triangle
//     mesh.vertices = (float *)std::memcpy(pool.allocate(render_mesh.vertices.size() * sizeof(float)), render_mesh.vertices.data(), render_mesh.vertices.size() * sizeof(float));
//     mesh.texcoords = (float *)std::memcpy(pool.allocate(render_mesh.texcoords.size() * sizeof(float)), render_mesh.texcoords.data(), render_mesh.texcoords.size() * sizeof(float));
//     mesh.normals = (float *)std::memcpy(pool.allocate(render_mesh.normals.size() * sizeof(float)), render_mesh.normals.data(), render_mesh.normals.size() * sizeof(float));
//     mesh.indices = NULL;
//     UploadMesh(&mesh, false);
//     return mesh;
// }

void copy_tangents_as_vec_4(float* dest, const std::vector<float>& tangents) {
    for (unsigned int i = 0; i < tangents.size() / 3; i++) {
        dest[4*i+0] = tangents[3 * i + 0]; // x
        dest[4*i+1] = tangents[3 * i + 1]; // y
        dest[4*i+2] = tangents[3 * i + 2]; // z
        dest[4*i+3] = 0.0f;                // w
    }
}

MeshGPU::MeshGPU(const RenderMesh& mesh) {
    // Allocate resources
    vertices = (float*) calloc(sizeof(float), mesh.vertices.size());
    texcoords = (float*) calloc(sizeof(float), mesh.texcoords.size());
    normals = (float*) calloc(sizeof(float), mesh.normals.size());
    tangents = (float*) calloc(sizeof(float), 4 * mesh.tangents.size() / 3);

    // Copy the data
    nVertices = mesh.vertices.size() / 3;
    vertices = (float *) std::memcpy(vertices, mesh.vertices.data(), mesh.vertices.size()*sizeof(float));
    texcoords = (float *) std::memcpy(texcoords, mesh.texcoords.data(), mesh.texcoords.size()*sizeof(float));
    normals = (float *) std::memcpy(normals, mesh.normals.data(), mesh.normals.size()*sizeof(float));
    copy_tangents_as_vec_4(tangents, mesh.tangents);
    if (!vertices) std::cerr << "MeshGPU::MeshGPU: vertices null" << std::endl;
    if (!texcoords) std::cerr << "MeshGPU::MeshGPU: texcoord null" << std::endl;
    if (!normals) std::cerr << "MeshGPU::MeshGPU: normals null" << std::endl;
    if (!tangents) std::cerr << "MeshGPU::MeshGPU: tangents null" << std::endl;

    // Upload the mesh to GPU
    Mesh raymesh = {0};
    raymesh.vertexCount = mesh.vertices.size() / 3; // 3 coordinates per vertex
    raymesh.triangleCount = raymesh.vertexCount / 3; // 3 vertices per triangle
    raymesh.vertices = vertices;
    raymesh.texcoords = texcoords;
    raymesh.normals = normals;
    raymesh.tangents = tangents;
    UploadMesh(&raymesh, true);
    // VBO indices handled by raylib
    verticesVBO = raymesh.vboId[0];
    texcoordsVBO = raymesh.vboId[1];
    normalsVBO = raymesh.vboId[2];
    // vboId[3] are colors which we do not support
    tangentsVBO = raymesh.vboId[4];
    VAO = raymesh.vaoId;
    RL_FREE(raymesh.vboId);
}

void MeshGPU::updateData(const RenderMesh& mesh) {
    // Copy the data from the render mesh
    vertices = (float*) std::memcpy(vertices, mesh.vertices.data(), mesh.vertices.size()*sizeof(float));
    normals = (float*) std::memcpy(normals, mesh.normals.data(), mesh.normals.size()*sizeof(float));
    copy_tangents_as_vec_4(tangents, mesh.tangents);

    // Update the VBOs in the GPU
    rlUpdateVertexBuffer(verticesVBO, vertices, mesh.vertices.size() * sizeof(float), 0);
    rlUpdateVertexBuffer(normalsVBO, normals, mesh.normals.size() * sizeof(float), 0);
    rlUpdateVertexBuffer(tangentsVBO, tangents, mesh.tangents.size() * sizeof(float), 0);
}

MeshGPU::~MeshGPU() {
    // Deallocate CPU resources
    free(vertices);
    free(texcoords);
    free(normals);
    free(tangents);

    // Unload the mesh buffers (deallocate in GPU)
    rlUnloadVertexArray(VAO);
    rlUnloadVertexBuffer(verticesVBO);
    rlUnloadVertexBuffer(texcoordsVBO);
    rlUnloadVertexBuffer(normalsVBO);
    rlUnloadVertexBuffer(tangentsVBO);
}

/**
 * MeshGPU handles the CPU and GPU memory for us, we use raylib's Mesh as a vector to
 * use raylib functionality of drawing meshes.
 *
 * Calling this function does not have a real cost as we are just copying data from one structure
 * to another and no memory management is being done.
 */
Mesh MeshGPUtoRaymesh(const MeshGPU& mesh, MemoryPool& pool) {
    Mesh raymesh = {0};
    raymesh.vertexCount = mesh.nVertices;
    raymesh.triangleCount = raymesh.vertexCount / 3; // 3 vertices per triangle
    raymesh.vertices = mesh.vertices;
    raymesh.texcoords = mesh.texcoords;
    raymesh.normals = mesh.normals;
    raymesh.tangents = mesh.tangents;

    // Copy the VAOs
    raymesh.vaoId = mesh.VAO;
    raymesh.vboId = (unsigned int*) pool.allocate(sizeof(unsigned int) * 3);
    raymesh.vboId[0] = mesh.verticesVBO;
    raymesh.vboId[1] = mesh.texcoordsVBO;
    raymesh.vboId[2] = mesh.normalsVBO;
    raymesh.vboId[4] = mesh.tangentsVBO;

    return raymesh;
}

// void UnloadGPUMesh(const Mesh& mesh) {
//     // Unload rlgl VAO and allocated VBOs
//     rlUnloadVertexArray(mesh.vaoId);

//     #define MAX_MESH_VERTEX_BUFFERS 7
//     if (mesh.vboId != NULL) for (int i = 0; i < MAX_MESH_VERTEX_BUFFERS; i++) rlUnloadVertexBuffer(mesh.vboId[i]);
//     RL_FREE(mesh.vboId);
//     #undef MAX_MESH_VERTEX_BUFFERS
// }

/**
 * INFO Global data for our renderer that will be used throughout the runtime of the application.
 *
 * NOTE We can not hold this variables inside of MandosViewer class, as they depend on raylib and MandosViewer should
 * be renderer agnostic by design.
 * REVIEW if this section goes out of hand, we should wrap all this variables in a struct and only have a global pointer to
 * a dinamically allocated instance of this struct. Memory handling should be done by MandosViewer.
 */

static Camera3D camera;
#define SPHERE_SUBDIVISIONS 30
static Model sphere_model;
static Material base_material;
static Shader base_shader;
static Shader normals_shader;
static Shader diffuse_shader;
static Shader axis_shader;
static Texture2D texture;
static Texture2D normalMapTexture;
static MeshGPU* axis3D = nullptr;

Material createMaterialFromShader(Shader shader) {
    Material material = { 0 };
    material.maps = (MaterialMap *)calloc(12, sizeof(MaterialMap));

    // Using rlgl default shader
    material.shader = shader;

    // Using rlgl default texture (1x1 pixel, UNCOMPRESSED_R8G8B8A8, 1 mipmap)
    material.maps[MATERIAL_MAP_DIFFUSE].texture = Texture2D( rlGetTextureIdDefault(), 1, 1, 1, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8 );

    material.maps[MATERIAL_MAP_DIFFUSE].color = WHITE;    // Diffuse color
    material.maps[MATERIAL_MAP_SPECULAR].color = WHITE;   // Specular color

    return material;
}

MandosViewer::MandosViewer() {
    InitWindow(initialScreenWidth, initialScreenHeight, "Mandos");
    SetWindowState(FLAG_WINDOW_RESIZABLE);
    rlDisableBackfaceCulling();
    camera = create_camera();
    sphere_model = LoadModelFromMesh(GenMeshSphere(0.1, SPHERE_SUBDIVISIONS, SPHERE_SUBDIVISIONS));
    SetTargetFPS(200);

    // Load basic lighting shader
    base_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/lighting.fs");
    base_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(base_shader, "viewPos");
    base_material = createMaterialFromShader(base_shader);

    // Other shaders
    normals_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/normals.fs");
    diffuse_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/picker.fs");
    axis_shader = LoadShader("resources/shaders/axis.vs", "resources/shaders/axis.fs");

    // Texture loading
    texture = LoadTexture("resources/textures/slab.png");
    normalMapTexture = LoadTexture("resources/textures/NormalMap.png");
    // SetMaterialTexture(&base_material, MATERIAL_MAP_DIFFUSE, texture);
    SetMaterialTexture(&base_material, MATERIAL_MAP_NORMAL, normalMapTexture);

    // Axis 3D
    RenderMesh axis3Drender = RenderMesh("resources/obj/axis.obj");
    axis3D = new MeshGPU(axis3Drender);

    // Get some required shader locations
    sphere_model.materialCount = 1;
    sphere_model.materials[0] = base_material;

}

MandosViewer::~MandosViewer() {
    UnloadModel(sphere_model); // Unloads also the material maps
    UnloadShader(base_shader);
    UnloadShader(normals_shader);
    UnloadShader(diffuse_shader);
    UnloadShader(axis_shader);
    UnloadTexture(texture);
    UnloadTexture(normalMapTexture);
    delete axis3D;

    CloseWindow(); // Destroys the opengl context
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

void draw_mesh_color(const Mat4& transform, const MeshGPU& mesh, MemoryPool& mem_pool, Color color) {
    base_material.maps[MATERIAL_MAP_DIFFUSE].color = color;
    Mesh raymesh = MeshGPUtoRaymesh(mesh, mem_pool);
    DrawMesh(raymesh, base_material, matrix_eigen_to_raylib(transform));
    base_material.maps[MATERIAL_MAP_DIFFUSE].color = WHITE;
}

void getPitchYaw(const Camera3D& camera, float& pitch, float& yaw) {
    // Calculate forward, right, and up vectors
    Vector3 forward = Vector3Subtract(camera.target, camera.position);
    forward = Vector3Normalize(forward);

    Vector3 right = Vector3CrossProduct(camera.up, forward);
    right = Vector3Normalize(right);

    // Calculate pitch
    pitch = asin(-forward.y);

    // Calculate yaw
    yaw = atan2(right.x, right.z);
}

void DrawAxis3D(MemoryPool& mem_pool) {
    Vector3 axisPosition = Vector3(30.0f, 15.0f, 0.0f);
    Camera3D axisCam = create_camera();
    axisCam.position = Vector3(0.0f, 0.0f, 10.0f);
    axisCam.projection = CAMERA_ORTHOGRAPHIC;
    Mat3 transform = Mat3::Identity()*2.f;
    Scalar pitch, yaw;
    getPitchYaw(camera, pitch, yaw);
    const Eigen::AngleAxis<Scalar> pitchRotation(pitch, Vec3(1.0f,0.0f,0.0f));
    const Eigen::AngleAxis<Scalar> yawRotation(-yaw - M_PI_2, Vec3(camera.up.x,camera.up.y,camera.up.z));
    transform = pitchRotation.toRotationMatrix() * yawRotation.toRotationMatrix() * transform;
    Mesh raymesh = MeshGPUtoRaymesh(*axis3D, mem_pool);

    Matrix rayTransform = {
    transform(0,0), transform(0,1), transform(0, 2), axisPosition.x,
    transform(1,0), transform(1,1), transform(1, 2), axisPosition.y,
    transform(2,0), transform(2,1), transform(2, 2), axisPosition.z,
    0.0f ,0.0f ,0.0f , 1.0f,
    };
    BeginMode3D(axisCam);
    base_material.shader = axis_shader;
    DrawMesh(raymesh, base_material, rayTransform);
    base_material.shader = base_shader;
    EndMode3D();
}

void MandosViewer::end_3D_mode() {
    EndMode3D();
    DrawAxis3D(mem_pool);
}

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
    DrawModel(sphere_model, vector3_eigen_to_raylib(position), 1, PARTICLE_COLOR);
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

void MandosViewer::draw_particle_indices(const Simulation& simulation, const PhysicsState& state) {
    for (size_t i = 0; i < simulation.simulables.particles.size(); i++) {
        const Particle particle = simulation.simulables.particles[i];
        const Vector3 position = vector3_eigen_to_raylib(particle.get_position(state.x));

        const unsigned int index = particle.index / 3;
        const Matrix matView = GetCameraMatrix(camera);
        const Vector3 particleCamPos = Vector3Transform(position, matView);

        // Filter out tags that are too far away or behind the camera
        const float thresholdDistance = 5.0f;
        if (particleCamPos.z > 0.0f or particleCamPos.z < -thresholdDistance) continue;

        Vector2 particleScreenSpacePosition = GetWorldToScreen(position, camera);
        DrawText(std::to_string(index).c_str(), (int)particleScreenSpacePosition.x - MeasureText(std::to_string(index).c_str(), 20)/2, (int)particleScreenSpacePosition.y, 20, BLACK);
    }
}
