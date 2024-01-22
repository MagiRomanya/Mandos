#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <stdlib.h>
#include <Eigen/Core>
#include <vector>

#include "fem_unit.hpp"
#include "imgui.h"
#include "linear_algebra.hpp"
#include "mandos.hpp"
#include "memory_pool.hpp"
#include "mesh.hpp"
#include "raylib.h"
#include "raylib_imgui.hpp"
#include "raymath.h"
#include "rlgl.h"
#include "rcamera.h"
#include "spring.hpp"
#include "viewmandos.hpp"
#include "utility_functions.hpp"

Color RB_COLOR = YELLOW;
Color FEM_COLOR = PINK;
Color PARTICLE_COLOR = BLUE;
Color MASS_SPRING_COLOR  = GREEN;
Color FROZEN_PARTICLE_COLOR = WHITE;

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
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type
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

/**
 * INFO Global data for our renderer that will be used throughout the runtime of the application.
 *
 * NOTE We can not hold this variables inside of MandosViewer class, as they depend on raylib and MandosViewer should
 * be renderer agnostic by design.
 */

#define SPHERE_SUBDIVISIONS 30
struct RenderState {
    Camera3D camera;
    Model sphere_model;
    Material base_material;

    Shader base_shader;
    Shader normals_shader;
    Shader pbr_shader;
    Shader axis_shader;

    Texture2D diffuseTexture;
    Texture2D normalMapTexture;

    MeshGPU* axis3D = nullptr;
};

static RenderState* renderState = nullptr;

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
    ImGuiInitialize();
    rlDisableBackfaceCulling();

    // Set up the render state
    renderState = new RenderState;
    renderState->camera = create_camera();
    renderState->sphere_model = LoadModelFromMesh(GenMeshSphere(0.1, SPHERE_SUBDIVISIONS, SPHERE_SUBDIVISIONS));
    SetTargetFPS(200);

    // Load basic lighting shader
    renderState->base_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/lighting.fs");
    renderState->base_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(renderState->base_shader, "viewPos");
    renderState->base_material = createMaterialFromShader(renderState->base_shader);

    // Other shaders
    renderState->normals_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/normals.fs");
    renderState->pbr_shader = LoadShader("resources/shaders/lighting.vs", "resources/shaders/pbr.fs");
    renderState->pbr_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(renderState->base_shader, "viewPos");
    renderState->axis_shader = LoadShader("resources/shaders/axis.vs", "resources/shaders/axis.fs");

    // Texture loading
    renderState->diffuseTexture = LoadTexture("resources/textures/terracota.png");
    renderState->normalMapTexture = LoadTexture("resources/textures/NormalMap.png");
    // SetMaterialTexture(&renderState->base_material, MATERIAL_MAP_DIFFUSE, renderState->diffuseTexture);
    // SetMaterialTexture(&renderState->base_material, MATERIAL_MAP_NORMAL, renderState->normalMapTexture);

    // Axis 3D
    RenderMesh axis3Drender = RenderMesh("resources/obj/axis.obj");
    renderState->axis3D = new MeshGPU(axis3Drender);

    // Get some required shader locations
    renderState->sphere_model.materialCount = 1;
    renderState->sphere_model.materials[0] = renderState->base_material;

}

MandosViewer::~MandosViewer() {
    UnloadModel(renderState->sphere_model); // Unloads also the material maps
    UnloadShader(renderState->base_shader);
    UnloadShader(renderState->normals_shader);
    UnloadShader(renderState->pbr_shader);
    UnloadShader(renderState->axis_shader);
    UnloadTexture(renderState->diffuseTexture);
    UnloadTexture(renderState->normalMapTexture);
    delete renderState->axis3D;

    ImGuiDeinitialize();
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
    ImGuiBeginDrawing();
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("Render")) {
            // Add menu items for "File"
            if (ImGui::MenuItem("New")) {
                // Handle "New" action
            }
            if (ImGui::MenuItem("Open")) {
                // Handle "Open" action
            }
            if (ImGui::MenuItem("Save")) {
                // Handle "Save" action
            }
            ImGui::EndMenu();
        }

        // Add more menus as needed...

        ImGui::EndMainMenuBar();
    }
    ImGuiEndDrawing();
    DrawFPS(GetScreenWidth()*0.95, GetScreenHeight()*0.05);
    EndDrawing();
}

void MandosViewer::begin_3D_mode() {
    BeginMode3D(renderState->camera);
    DrawGrid(30, 1.0f);
}

void draw_mesh_color(const Mat4& transform, const MeshGPU& mesh, MemoryPool& mem_pool, Color color) {
    renderState->base_material.maps[MATERIAL_MAP_DIFFUSE].color = color;
    Mesh raymesh = MeshGPUtoRaymesh(mesh, mem_pool);
    DrawMesh(raymesh, renderState->base_material, matrix_eigen_to_raylib(transform));
    renderState->base_material.maps[MATERIAL_MAP_DIFFUSE].color = WHITE;
}

void getPitchYaw(const Camera3D& camera, float& pitch, float& yaw) {
    // Calculate forward, right, and up vectors
    Vector3 forward = Vector3Subtract(camera.target, camera.position);
    forward = Vector3Normalize(forward);

    Vector3 right = Vector3CrossProduct(camera.up, forward);
    right = Vector3Normalize(right);

    pitch = asin(-forward.y);
    yaw = atan2(right.x, right.z);
}

void DrawAxis3D(MemoryPool& mem_pool) {
    Vector3 axisPosition = Vector3(30.0f, 15.0f, 0.0f);
    Camera3D axisCam = create_camera();
    axisCam.position = Vector3(0.0f, 0.0f, 10.0f);
    axisCam.projection = CAMERA_ORTHOGRAPHIC;
    Mat3 transform = Mat3::Identity()*2.f;
    Scalar pitch, yaw;
    getPitchYaw(renderState->camera, pitch, yaw);
    const Eigen::AngleAxis<Scalar> pitchRotation(pitch, Vec3(1.0f,0.0f,0.0f));
    const Eigen::AngleAxis<Scalar> yawRotation(-yaw - M_PI_2, Vec3(renderState->camera.up.x,renderState->camera.up.y,renderState->camera.up.z));
    transform = pitchRotation.toRotationMatrix() * yawRotation.toRotationMatrix() * transform;
    Mesh raymesh = MeshGPUtoRaymesh(*renderState->axis3D, mem_pool);

    Matrix rayTransform = {
    transform(0,0), transform(0,1), transform(0, 2), axisPosition.x,
    transform(1,0), transform(1,1), transform(1, 2), axisPosition.y,
    transform(2,0), transform(2,1), transform(2, 2), axisPosition.z,
    0.0f ,0.0f ,0.0f , 1.0f,
    };
    BeginMode3D(axisCam);
    renderState->base_material.shader = renderState->axis_shader;
    DrawMesh(raymesh, renderState->base_material, rayTransform);
    renderState->base_material.shader = renderState->base_shader;
    EndMode3D();
}

void MandosViewer::end_3D_mode() {
    EndMode3D();
    DrawAxis3D(mem_pool);
}

void myUpdateCamera(Camera3D& camera) {
    const float CAMERA_SENSITIVITY = 0.01f;
    const float CAMERA_MOVE_SPEED = 0.01f;
    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
        const Vector2 delta = GetMouseDelta();
        if (IsKeyDown(KEY_LEFT_SHIFT)) {
            CameraMoveRight(&camera, -CAMERA_MOVE_SPEED*delta.x, true);
            CameraMoveUp(&camera, CAMERA_MOVE_SPEED*delta.y);
        }
        else {
            CameraYaw(&camera, -delta.x*CAMERA_SENSITIVITY, true);
            CameraPitch(&camera, -delta.y*CAMERA_SENSITIVITY, true, true, false);
        }
    }
    // Zoom target distance
    CameraMoveToTarget(&camera, -GetMouseWheelMove());
}

void MandosViewer::update_camera() {
    myUpdateCamera(renderState->camera);

    // Update the shaders with the camera view vector (points towards { 0.0f, 0.0f, 0.0f })
    float cameraPos[3] = { renderState->camera.position.x, renderState->camera.position.y, renderState->camera.position.z };
    SetShaderValue(renderState->base_shader, renderState->base_shader.locs[SHADER_LOC_VECTOR_VIEW], cameraPos, SHADER_UNIFORM_VEC3);
    SetShaderValue(renderState->pbr_shader, renderState->pbr_shader.locs[SHADER_LOC_VECTOR_VIEW], cameraPos, SHADER_UNIFORM_VEC3);
}

bool MandosViewer::is_key_pressed(int Key) { return IsKeyPressed(Key); }

void MandosViewer::draw_particle(const ParticleHandle& particle, const PhysicsState& state) {
    const Vec3 position = particle.particle.get_position(state.x);
    DrawModel(renderState->sphere_model, vector3_eigen_to_raylib(position), 1, PARTICLE_COLOR);
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
        DrawModel(renderState->sphere_model, vector3_eigen_to_raylib(position), 1, color);
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
        const Matrix matView = GetCameraMatrix(renderState->camera);
        const Vector3 particleCamPos = Vector3Transform(position, matView);

        // Filter out tags that are too far away or behind the camera
        const float thresholdDistance = 5.0f;
        if (particleCamPos.z > 0.0f or particleCamPos.z < -thresholdDistance) continue;

        Vector2 particleScreenSpacePosition = GetWorldToScreen(position, renderState->camera);
        DrawText(std::to_string(index).c_str(), (int)particleScreenSpacePosition.x - MeasureText(std::to_string(index).c_str(), 20)/2, (int)particleScreenSpacePosition.y, 20, BLACK);
    }
}
