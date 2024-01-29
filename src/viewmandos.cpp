#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <stdlib.h>
#include <Eigen/Core>
#include <vector>
#include <imgui.h>
#include <ImGuizmo.h>

#include "fem_unit.hpp"
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
#include "gl_extensions.hpp"

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

/**
 * INFO Global data for our renderer that will be used throughout the runtime of the application.
 *
 * NOTE We can not hold this variables inside of MandosViewer class, as they depend on raylib and MandosViewer should
 * be renderer agnostic by design.
 */

#define SPHERE_SUBDIVISIONS 30
#define SHADER_LIST \
    SHADER(normals_shader, "resources/shaders/lighting.vs", "resources/shaders/normals.fs") \
    SHADER(pbr_shader, "resources/shaders/lighting.vs", "resources/shaders/pbr.fs") \
    SHADER(bling_phong_shader, "resources/shaders/lighting.vs", "resources/shaders/lighting.fs") \
    SHADER(bling_phong_texture_shader, "resources/shaders/lighting.vs", "resources/shaders/bling-phong-textures.fs") \
    SHADER(axis_shader, "resources/shaders/axis.vs", "resources/shaders/axis.fs") \
    SHADER(instancing_shader, "resources/shaders/lighting_instancing.vs", "resources/shaders/lighting.fs") \
    SHADER(solid_shader, "resources/shaders/lighting.vs", "resources/shaders/diffuse_solid_color.fs")

struct RenderState {
    void initialize();
    void deinitialize();

    // MISC
    // -------------------------------------------------------------
    Camera3D camera;
    Material base_material;

    // SHADERS
    // -------------------------------------------------------------
    Shader base_shader; // DO NOT INITIALIZE OR DESTROY
#define SHADER(var, vspath, fspath) Shader var;
    SHADER_LIST
#undef SHADER

    // TEXTURES
    // -------------------------------------------------------------
    Texture2D diffuseTexture;
    Texture2D normalMapTexture;
    RenderTexture2D screenFBO;

    // GEOMETRY
    // -------------------------------------------------------------
    Mesh sphere_mesh;
    MeshGPU* axis3D = nullptr;
    LinesGPU* lines = nullptr;

    // GUI
    // -------------------------------------------------------------
    bool renderFrameOpen = false;
    bool cameraFrameOpen = false;
    bool ImGuiLogsOpen = false;
    bool RaylibLogsOpen = false;
    std::vector<std::pair<int, std::string>> raylibLogs;
};

#define SHADER_LOC_SLICE_PLANE SHADER_LOC_COLOR_AMBIENT

void RenderState::initialize() {
    camera = create_camera();

    // Load shaders
#define SHADER(var, vspath, fspath) var = LoadShader(vspath, fspath);
    SHADER_LIST
#undef SHADER

    // Shaders set uniforms
    bling_phong_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(bling_phong_shader, "viewPos");
    bling_phong_shader.locs[SHADER_LOC_SLICE_PLANE] = GetShaderLocation(bling_phong_shader, "slicePlane");
    bling_phong_texture_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(bling_phong_texture_shader, "viewPos");
    bling_phong_texture_shader.locs[SHADER_LOC_SLICE_PLANE] = GetShaderLocation(bling_phong_texture_shader, "slicePlane");
    pbr_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(pbr_shader, "viewPos");
    pbr_shader.locs[SHADER_LOC_SLICE_PLANE] = GetShaderLocation(pbr_shader, "slicePlane");
    instancing_shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(instancing_shader, "viewPos");
    instancing_shader.locs[SHADER_LOC_MATRIX_MODEL] = GetShaderLocationAttrib(instancing_shader, "instanceTransform");
    solid_shader.locs[SHADER_LOC_SLICE_PLANE] = GetShaderLocation(solid_shader, "slicePlane");
    normals_shader.locs[SHADER_LOC_SLICE_PLANE] = GetShaderLocation(normals_shader, "slicePlane");

    // Set sefault shader
    base_shader = bling_phong_shader;
    base_material = createMaterialFromShader(base_shader);

    // Texture loading
    diffuseTexture = LoadTexture("resources/textures/terracota.png");
    normalMapTexture = LoadTexture("resources/textures/NormalMap.png");
    SetMaterialTexture(&base_material, MATERIAL_MAP_DIFFUSE, diffuseTexture);
    SetMaterialTexture(&base_material, MATERIAL_MAP_NORMAL, normalMapTexture);

    // Meshes
    sphere_mesh = GenMeshSphere(0.1, SPHERE_SUBDIVISIONS, SPHERE_SUBDIVISIONS);
    RenderMesh axis3Drender = RenderMesh("resources/obj/axis.obj");
    axis3D = new MeshGPU(axis3Drender);
    std::vector<float> vertices = {0.0f, 1.0f, 4.0f, 2.0f};
    lines = new LinesGPU(vertices);

    // Set up the FBO
    screenFBO = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());
}

void RenderState::deinitialize() {
    // Unload geometry
    UnloadMesh(sphere_mesh);
    delete axis3D;
    delete lines;

    // Unload shaders
#define SHADER(var, vspath, fspath) UnloadShader(var);
    SHADER_LIST
#undef SHADER

    // Unload textures
    UnloadTexture(diffuseTexture);
    UnloadTexture(normalMapTexture);
    UnloadRenderTexture(screenFBO);
}

static RenderState* renderState = nullptr;

void raylib_to_imgui_tracelog_callback(int logLevel, const char *text, va_list args) {
    std::string logString;
    switch (logLevel)
    {
        case LOG_TRACE: logString = "TRACE: "; break;
        case LOG_DEBUG: logString = "DEBUG: "; break;
        case LOG_INFO: logString = "INFO: "; break;
        case LOG_WARNING: logString = "WARNING: "; break;
        case LOG_ERROR: logString = "ERROR: "; break;
        case LOG_FATAL: logString = "FATAL: "; break;
        default: break;
    }
    const int bufferSize = 512;
    char buffer[bufferSize];
    vsnprintf(buffer, bufferSize, text, args);
    logString += std::string(buffer);
    logString += "\n";

    renderState->raylibLogs.emplace_back(logLevel, logString);
    std::cout << logString << std::flush;
}

MandosViewer::MandosViewer() {
    InitWindow(initialScreenWidth, initialScreenHeight, "Mandos");
    ImGuiInitialize();
    SetTraceLogCallback(raylib_to_imgui_tracelog_callback);
    SetWindowState(FLAG_WINDOW_RESIZABLE);
    rlDisableBackfaceCulling();

    // Set up the render state
    renderState = new RenderState;
    renderState->initialize();

    SetTargetFPS(200);
}

MandosViewer::~MandosViewer() {
    SetTraceLogCallback(NULL);
    renderState->deinitialize();
    delete renderState;
    ImGuiDeinitialize();
    CloseWindow(); // Destroys the opengl context
}

bool MandosViewer::window_should_close() {
    return WindowShouldClose();
}

void MandosViewer::begin_drawing() {
    mem_pool.reset();
    BeginDrawing();
    ClearBackground(WHITE);
    BeginTextureMode(renderState->screenFBO);
    // Update slice plane uniform
    SetShaderValue(renderState->base_shader, renderState->base_shader.locs[SHADER_LOC_SLICE_PLANE], slicePlane.data(), SHADER_UNIFORM_VEC4);
    SetShaderValue(renderState->solid_shader, renderState->solid_shader.locs[SHADER_LOC_SLICE_PLANE], slicePlane.data(), SHADER_UNIFORM_VEC4);
    ClearBackground(RAYWHITE);
}

inline void getFloatsFromColor(float* fc, Color color) {
    fc[0] = (float) color.r / 255.f;
    fc[1] = (float) color.g / 255.f;
    fc[2] = (float) color.b / 255.f;
    fc[3] = (float) color.a / 255.f;
}

inline Color getColorFromFloats(float* fc) {
    return Color((unsigned char) 255.f * fc[0], (unsigned char) 255.f * fc[1], (unsigned char) 255.f * fc[2], (unsigned char) 255.f * fc[3]);
}

inline void GUI_control_color(Color& color, std::string name) {
    float floatColor[4];
    getFloatsFromColor(floatColor, color);
    if (ImGui::ColorEdit4(name.c_str(), floatColor)) {
        color = getColorFromFloats(floatColor);
    }
}

inline void GUI_ShowRaylibLogs() {
    enum myTraceLevel {
    myLOG_NONE = 0,
    myLOG_WARNING = 1 << 1,
    myLOG_INFO = 1 << 2,
    myLOG_DEBUG = 1 << 3,
    myLOG_ERROR = 1 << 4,
    myLOG_ALL = myLOG_WARNING | myLOG_INFO | myLOG_DEBUG | myLOG_ERROR,
    };
    static int traceLevel = 0;
    ImGui::CheckboxFlags("All", &traceLevel, myLOG_ALL);
    ImGui::SameLine(); ImGui::CheckboxFlags("Warning", &traceLevel, myLOG_WARNING);
    ImGui::SameLine(); ImGui::CheckboxFlags("Info",  &traceLevel, myLOG_INFO);
    ImGui::SameLine(); ImGui::CheckboxFlags("Debug", &traceLevel, myLOG_DEBUG);
    ImGui::SameLine(); ImGui::CheckboxFlags("Error", &traceLevel, myLOG_ERROR);

    ImGui::BeginChild("##log", ImVec2(0.0f, 0.0f), ImGuiChildFlags_Border, ImGuiWindowFlags_AlwaysVerticalScrollbar | ImGuiWindowFlags_AlwaysHorizontalScrollbar);
    ImGuiListClipper clipper;
    clipper.Begin(renderState->raylibLogs.size(), ImGui::GetTextLineHeightWithSpacing());
    while (clipper.Step()) {
        for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; i++) {

            myTraceLevel currentLevel = (myTraceLevel)traceLevel;
            if (currentLevel == myLOG_ALL) {
                ImGui::Text("%s", renderState->raylibLogs[i].second.c_str());
            }
            else {
                if ((currentLevel & myLOG_WARNING) && renderState->raylibLogs[i].first == LOG_WARNING) {
                    ImGui::Text("%s", renderState->raylibLogs[i].second.c_str());
                }
                if ((currentLevel & myLOG_INFO) && renderState->raylibLogs[i].first == LOG_INFO) {
                    ImGui::Text("%s", renderState->raylibLogs[i].second.c_str());
                }

                if ((currentLevel & myLOG_DEBUG) && renderState->raylibLogs[i].first == LOG_DEBUG) {
                    ImGui::Text("%s", renderState->raylibLogs[i].second.c_str());
                }
                if ((currentLevel & myLOG_ERROR) && renderState->raylibLogs[i].first == LOG_ERROR) {
                    ImGui::Text("%s", renderState->raylibLogs[i].second.c_str());
                }

            }
        }
    }
    if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
        ImGui::SetScrollHereY(1.0f);
    ImGui::EndChild();
}

static inline void drawSlicingPlaneGuizmo(Vec4& slicePlaneVec, ImGuizmo::OPERATION operation) {
    // PLANE GUIZMO
    ImGuizmo::SetOrthographic(false);
    ImGuizmo::SetDrawlist();
    const float screenWidth = ImGui::GetContentRegionAvail().x;
    const float screenHeight = ImGui::GetContentRegionAvail().y;
    const float aspectRatio = screenWidth / screenHeight;
    const ImVec2 windowPos = ImGui::GetWindowPos();

    // set the screen dimensions
    ImGuizmo::SetRect(windowPos.x, windowPos.y, screenWidth, screenHeight);

    // We need to transpose the matrices because of the Matrix struct memory order of the components (column-major)
    const Matrix matView = MatrixTranspose(GetCameraViewMatrix(&renderState->camera));
    const Matrix matProj = MatrixTranspose(GetCameraProjectionMatrix(&renderState->camera, aspectRatio));
    static Matrix planeTransform = MatrixIdentity();
    ImGuizmo::Manipulate((float *)&matView, (float *)&matProj, operation, ImGuizmo::LOCAL, (float *)&planeTransform);
    // Keep in mind Transform is transposed!
    Mat3 planeRotation;
    planeRotation <<
        planeTransform.m0, planeTransform.m1, planeTransform.m2,
        planeTransform.m4, planeTransform.m5, planeTransform.m6,
        planeTransform.m8, planeTransform.m9, planeTransform.m10;
    Vec3 planeTranslation = Vec3(planeTransform.m3, planeTransform.m7, planeTransform.m11);

    const Vec3 initialPlaneNormal = Vec3(0.0f,1.0f,0.0f);
    const Vec3 planeNormal = planeRotation * initialPlaneNormal;
    const float D = - planeTranslation.dot(planeNormal);
    slicePlaneVec = Vec4(planeNormal.x(), planeNormal.y(), planeNormal.z(), D);
}

void myUpdateCamera(Camera3D& camera) {
    // if (ImGui::GetIO().WantCaptureMouse) return;
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

    // Update the shaders with the camera view vector (points towards { 0.0f, 0.0f, 0.0f })
    float cameraPos[3] = { renderState->camera.position.x, renderState->camera.position.y, renderState->camera.position.z };
    SetShaderValue(renderState->base_shader, renderState->base_shader.locs[SHADER_LOC_VECTOR_VIEW], cameraPos, SHADER_UNIFORM_VEC3);
    SetShaderValue(renderState->instancing_shader, renderState->instancing_shader.locs[SHADER_LOC_VECTOR_VIEW], cameraPos, SHADER_UNIFORM_VEC3);
}

void MandosViewer::drawSimulationVisualizationWindow() {
    bool VisualizationOpen = true;
    ImGuiWindowFlags flags = ImGuiWindowFlags_NoCollapse
        | ImGuiWindowFlags_NoScrollWithMouse
        | ImGuiWindowFlags_NoScrollbar
        ;
    ImGui::Begin("Simulation view", &VisualizationOpen, flags);
    if (ImGui::IsWindowHovered() && not ImGuizmo::IsOver()) {
        myUpdateCamera(renderState->camera);
    }
    const float window_width = ImGui::GetContentRegionAvail().x;
    const float window_height = ImGui::GetContentRegionAvail().y;
    simulationViewerWidth = window_width;
    simulationViewerHeight = window_height;

    // Resizing the FBO
    UpdateRenderTexture2D(renderState->screenFBO, window_width, window_height);
    rlViewport(0, 0, window_width, window_height);

    // Screen position of the window
    ImVec2 pos = ImGui::GetCursorScreenPos();

    // Add the Simulation Render to the window
    ImGui::GetWindowDrawList()->AddImage(
        reinterpret_cast<ImTextureID>(renderState->screenFBO.texture.id),
        ImVec2(pos.x, pos.y),
        ImVec2(pos.x + window_width, pos.y + window_height),
        ImVec2(0, 1),
        ImVec2(1, 0)
        );

    if (enable_slice_plane) {
        if (!slicePlaneGuizmoToggle)
            drawSlicingPlaneGuizmo(slicePlane, ImGuizmo::OPERATION::ROTATE);
        else
            drawSlicingPlaneGuizmo(slicePlane, ImGuizmo::OPERATION::TRANSLATE);
    }
    ImGui::End();
}

void MandosViewer::drawGUI() {
    ImGuiBeginDrawing();
    drawSimulationVisualizationWindow();

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("Menu")) {
            if (ImGui::MenuItem("Render")) {
                renderState->renderFrameOpen = true;
            }
            if (ImGui::MenuItem("Camera")) {
                renderState->cameraFrameOpen = true;
            }
            if (ImGui::BeginMenu("Logs")) {
                if (ImGui::MenuItem("Raylib Logs"))
                    renderState->RaylibLogsOpen = true;
                if (ImGui::MenuItem("ImGui Logs"))
                    renderState->ImGuiLogsOpen = true;
                ImGui::EndMenu();
            }
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
    if (renderState->renderFrameOpen) {
        ImGui::Begin("Render", &renderState->renderFrameOpen);

        ImGui::SeparatorText("Simulable colors");
        GUI_control_color(PARTICLE_COLOR, "Particle color");
        GUI_control_color(RB_COLOR, "Rigid Body color");
        GUI_control_color(MASS_SPRING_COLOR, "Mass spring color");
        GUI_control_color(FEM_COLOR, "FEM color");

        ImGui::SeparatorText("Simulable shaders");
        ImVec2 buttonBox = ImVec2(ImGui::GetContentRegionAvail().x, 30);
        if (ImGui::Button("Bling Phong shader", buttonBox)) {
            renderState->base_shader = renderState->bling_phong_shader;
        }
        if (ImGui::Button("Bling Phong texture shader", buttonBox)) {
            renderState->base_shader = renderState->bling_phong_texture_shader;
        }
        if (ImGui::Button("PBR", buttonBox)) {
            renderState->base_shader = renderState->pbr_shader;
        }
        if (ImGui::Button("Normals shader", buttonBox)) {
            renderState->base_shader = renderState->normals_shader;
        }
        if (ImGui::Button("Solid Color", buttonBox)) {
            renderState->base_shader = renderState->solid_shader;
        }
        ImGui::SeparatorText("Simulation state visualization");
        ImGui::Checkbox("Render simulable meshes", &enable_draw_simulable_meshes);
        ImGui::Checkbox("Render particles", &enable_draw_particles);
        ImGui::Checkbox("Render springs", &enable_draw_springs);
        ImGui::Checkbox("Render FEM tetrahedrons", &enable_draw_fem_tetrahedrons);
        ImGui::Checkbox("Render particle indices", &enable_draw_particle_indices);

        ImGui::SeparatorText("Slice plane");
        ImGui::Checkbox("Enable Slice plane", &enable_slice_plane);
        if (enable_slice_plane) {
            EnableUserDefinedClipping();
            if (ImGui::Button("Toggle Translation Rotation", buttonBox)) {
                slicePlaneGuizmoToggle = not slicePlaneGuizmoToggle;
            }
            ImGui::InputFloat4("Slice plane", slicePlane.data());
        }
        else {
            DisableUserDefinedClipping();
        }

        ImGui::End();
    }
    if (renderState->ImGuiLogsOpen) {
        ImGui::ShowDebugLogWindow(&renderState->ImGuiLogsOpen);
    }
    if (renderState->RaylibLogsOpen) {
        ImGui::Begin("Raylib Logs", &renderState->RaylibLogsOpen);
        GUI_ShowRaylibLogs();
        ImGui::End();
    }
    ImGuiEndDrawing();
}

void MandosViewer::end_drawing() {
    EndTextureMode();
    drawGUI();
    EndDrawing();
}

void MandosViewer::begin_3D_mode() {
    BeginMode3D(renderState->camera);
    DrawGrid(100, 1.0f);
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
    if (enable_draw_particle_indices && SavedSim && SavedState) {
        draw_particle_indices(*SavedSim, *SavedState);
    }
}

bool MandosViewer::is_key_pressed(int Key) { return IsKeyPressed(Key); }

void MandosViewer::draw_particle(const ParticleHandle& particle, const PhysicsState& state) {
    if (!enable_draw_simulable_meshes) return;
    const Vec3 position = particle.particle.get_position(state.x);
    Matrix transform = MatrixTranslate(position.x(), position.y(), position.z());
    DrawMesh(renderState->sphere_mesh, renderState->base_material, transform);
}

void MandosViewer::draw_rigid_body(const RigidBodyHandle& rb, const PhysicsState& state, const MeshGPU& mesh) {
    if (!enable_draw_simulable_meshes) return;
    draw_mesh_color(rb.get_transformation_matrix(state), mesh, mem_pool, RB_COLOR);
}
void MandosViewer::draw_mesh(const Mat4& transform, const MeshGPU& mesh) {
    if (!enable_draw_simulable_meshes) return;
    draw_mesh_color(transform, mesh, mem_pool, WHITE);
}

void MandosViewer::draw_springs(const Simulation& simulation, const PhysicsState& state) {
    std::vector<float> vertices;
    for (size_t i = 0; i < simulation.energies.particle_springs.size(); i++) {
        ParticleSpring s = simulation.energies.particle_springs[i];
        Vec3 x1 = s.p1.get_position(state.x);
        Vec3 x2 = s.p2.get_position(state.x);
        vertices.push_back(x1.x());
        vertices.push_back(x1.y());
        vertices.push_back(x1.z());
        vertices.push_back(x2.x());
        vertices.push_back(x2.y());
        vertices.push_back(x2.z());
    }

    Material material = createMaterialFromShader(renderState->solid_shader);
    material.maps[MATERIAL_MAP_DIFFUSE].color = RED;
    renderState->lines->drawLines(material, vertices);
    free(material.maps);
}

void MandosViewer::draw_particles(const Simulation& simulation, const PhysicsState& state) {
    std::vector<Matrix> transforms;
    transforms.reserve(simulation.simulables.particles.size());

    for (size_t i = 0; i < simulation.simulables.particles.size(); i++) {
        Particle particle = simulation.simulables.particles[i];
        Vec3 position = particle.get_position(state.x);
        Matrix transform = MatrixTranslate(position.x(), position.y(), position.z());
        transforms.push_back(transform);
    }
    Material matInstances = LoadMaterialDefault();
    matInstances.shader = renderState->instancing_shader;
    matInstances.maps[MATERIAL_MAP_DIFFUSE].color = PARTICLE_COLOR;
    DrawMeshInstanced(renderState->sphere_mesh, matInstances, transforms.data(), transforms.size());
}

void MandosViewer::draw_FEM(const FEMHandle& fem, const PhysicsState& state, MeshGPU& gpuMesh, RenderMesh& renderMesh, SimulationMesh& simMesh) {
    if (!enable_draw_simulable_meshes) return;
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
    if (!enable_draw_simulable_meshes) return;
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
    std::vector<float> vertices;

#define MAT(type, name) \
    for (unsigned int i = 0; i < simulation.energies.fem_elements_##name.size(); i++) { \
        const FEM_Element3D<type>& e = simulation.energies.fem_elements_##name[i]; \
        const Vec3& x1 = e.p1.get_position(state.x);  \
        const Vec3& x2 = e.p2.get_position(state.x);  \
        const Vec3& x3 = e.p3.get_position(state.x);  \
        const Vec3& x4 = e.p4.get_position(state.x);  \
        vertices.insert(vertices.end(), {x1.x(), x1.y(), x1.z()}); \
        vertices.insert(vertices.end(), {x2.x(), x2.y(), x2.z()}); \
        vertices.insert(vertices.end(), {x1.x(), x1.y(), x1.z()}); \
        vertices.insert(vertices.end(), {x3.x(), x3.y(), x3.z()}); \
        vertices.insert(vertices.end(), {x1.x(), x1.y(), x1.z()}); \
        vertices.insert(vertices.end(), {x4.x(), x4.y(), x4.z()}); \
        vertices.insert(vertices.end(), {x2.x(), x2.y(), x2.z()}); \
        vertices.insert(vertices.end(), {x4.x(), x4.y(), x4.z()}); \
        vertices.insert(vertices.end(), {x2.x(), x2.y(), x2.z()}); \
        vertices.insert(vertices.end(), {x3.x(), x3.y(), x3.z()}); \
        vertices.insert(vertices.end(), {x3.x(), x3.y(), x3.z()}); \
        vertices.insert(vertices.end(), {x4.x(), x4.y(), x4.z()}); \
    }
    FEM_MATERIAL_MEMBERS
#undef MAT

    Material material = createMaterialFromShader(renderState->solid_shader);
    material.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    renderState->lines->drawLines(material, vertices);

    free(material.maps);
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

        Vector2 particleScreenSpacePosition = GetWorldToScreenEx(position, renderState->camera, simulationViewerWidth, simulationViewerHeight);
        DrawText(std::to_string(index).c_str(), (int)particleScreenSpacePosition.x - MeasureText(std::to_string(index).c_str(), 20)/2, (int)particleScreenSpacePosition.y, 20, BLACK);
    }
}

void MandosViewer::draw_simulation_state(const Simulation& simulation, const PhysicsState& state) {
    SavedSim = &simulation;
    SavedState = &state;
    if (enable_draw_springs)
        draw_springs(simulation, state);
    if (enable_draw_particles)
        draw_particles(simulation, state);
    if (enable_draw_fem_tetrahedrons)
        draw_FEM_tetrahedrons(simulation, state);
}
