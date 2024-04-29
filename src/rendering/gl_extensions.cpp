#include <cassert>
#include <iostream>
#include <vector>

#include "external/glad.h"
#include "gl_extensions.hpp"
#include "raylib.h"
#include "raymath.h"
#include "viewmandos.hpp"
#include <rlgl.h>

// OPENGL ERROR HANDLING
// ----------------------------------------------------------------
static void GLClearError() {
    while (glGetError() != GL_NO_ERROR);
}

static bool GLLogCall(std::string call_name, std::string file, int line) {
    while (GLenum error = glGetError()) {
        std::cerr << "[OpenGL ERROR] (" << error << ") in file "
            << file << ", and line " << line <<  "\n"
            << call_name << std::endl;
        return false;
    }
    return true;
}

#define GL_CALL(x) \
    GLClearError(); \
    x; \
    assert(GLLogCall(#x, __FILE__, __LINE__))
// ----------------------------------------------------------------


void EnableUserDefinedClipping() {
    GL_CALL(glEnable(GL_CLIP_DISTANCE0));
}

void DisableUserDefinedClipping() {
    GL_CALL(glDisable(GL_CLIP_DISTANCE0));
}


// Multisampled framebuffer, color and depth attachements
unsigned int msFBO, msRBO, msTexture;

int ANTIALIASING_SAMPLES = 4;

void InitializeMultisampleFramebuffer() {
    int max_samples;
    glGetIntegerv(GL_MAX_COLOR_TEXTURE_SAMPLES, &max_samples);
    if (max_samples < ANTIALIASING_SAMPLES) {
        std::cerr << "ERROR:FRAMEBUFFER:: The number of color texture samples is too big and not supported. Using the maximum amount of samples."  << std::endl;
        ANTIALIASING_SAMPLES = max_samples;
    }

    GL_CALL(glEnable(GL_MULTISAMPLE));

    // Create the multisapled framebuffer
    GL_CALL(glGenFramebuffers(1, &msFBO));
    GL_CALL(glBindFramebuffer(GL_FRAMEBUFFER, msFBO));

    // Create the multisampled color texture
    GL_CALL(glGenTextures(1, &msTexture));
    GL_CALL(glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, msTexture));
    GL_CALL(glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, ANTIALIASING_SAMPLES, GL_RGBA, GetScreenWidth(), GetScreenHeight(), 1));
    GL_CALL(glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0));
    GL_CALL(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, msTexture, 0));

    // Create the multisampled depthbuffer
    GL_CALL(glGenRenderbuffers(1, &msRBO));
    GL_CALL(glBindRenderbuffer(GL_RENDERBUFFER, msRBO));
    GL_CALL(glRenderbufferStorageMultisample(GL_RENDERBUFFER, ANTIALIASING_SAMPLES, GL_DEPTH24_STENCIL8, GetScreenWidth(), GetScreenHeight()));
    GL_CALL(glBindRenderbuffer(GL_RENDERBUFFER, 0));
    GL_CALL(glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, msRBO));

    // Check weather the framebuffer has the correct attachments
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cerr << "ERROR:FRAMEBUFFER:: Multisampled frame buffer is not complete!" << std::endl;
    }

    GL_CALL(glBindFramebuffer(GL_FRAMEBUFFER, 0));
}

void BindMultisampleFramebuffer() {
    GL_CALL(glBindFramebuffer(GL_FRAMEBUFFER, msFBO));
}

void BlitMultisampleFramebuffer(const RenderTexture destination) {
    const int width = destination.texture.width;
    const int height = destination.texture.height;
    GL_CALL(glBindFramebuffer(GL_READ_FRAMEBUFFER, msFBO));
    GL_CALL(glBindFramebuffer(GL_DRAW_FRAMEBUFFER, destination.id));
    GL_CALL(glBlitFramebuffer(0,0, width, height, 0,0, width, height, GL_COLOR_BUFFER_BIT, GL_NEAREST));
}

void UpdateRenderTexture2D(RenderTexture2D& fbo, int width, int height) {
    // Only execute the function if the width and height of the FBO are different from the new ones.
    if (width <= 0 or height <= 0) return;
    if (fbo.texture.width == width && fbo.texture.height == height) return;

    // Resize FBO
    GL_CALL(glBindFramebuffer(GL_FRAMEBUFFER, fbo.id));

    GL_CALL(glBindTexture(GL_TEXTURE_2D, fbo.texture.id));
    GL_CALL(glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL));
    GL_CALL(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
    GL_CALL(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
    GL_CALL(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo.texture.id, 0));

    // Framebuffer depth texture attachment (it is a renderbuffer)
    GL_CALL(glBindRenderbuffer(GL_RENDERBUFFER, fbo.depth.id));
    GL_CALL(glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height));
    GL_CALL(glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, fbo.depth.id));

    GL_CALL(glBindFramebuffer(GL_FRAMEBUFFER, 0));

    // Resize Multisample FBO for anti aliasing if it exists!
    GL_CALL(glBindFramebuffer(GL_FRAMEBUFFER, msFBO));
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE) {
        // Color attachment
        GL_CALL(glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, msTexture));
        GL_CALL(glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, ANTIALIASING_SAMPLES, GL_RGBA, width, height, GL_TRUE));
        GL_CALL(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, msTexture, 0));

        // Depth attachment
        GL_CALL(glBindRenderbuffer(GL_RENDERBUFFER, msRBO));
        GL_CALL(glRenderbufferStorageMultisample(GL_RENDERBUFFER, ANTIALIASING_SAMPLES, GL_DEPTH24_STENCIL8, width, height));
        GL_CALL(glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, msRBO));

    }
    GL_CALL(glBindFramebuffer(GL_FRAMEBUFFER, 0));

    // Update raylib FBO width and height
    rlSetFramebufferWidth(width);
    rlSetFramebufferHeight(height);
    fbo.texture.width = width;
    fbo.texture.height = height;
    fbo.depth.width = width;
    fbo.depth.height = height;
}

Texture1D CreateTexture1D() {
    Texture1D texture;
    GL_CALL(glGenTextures(1, &texture.id));
    return texture;
}

void UseTexture1D(Texture1D& texture) {
    GL_CALL(glActiveTexture(GL_TEXTURE0));
    GL_CALL(glBindTexture(GL_TEXTURE_1D, texture.id));
}

void UpdateTexture1D(Texture1D texture, const std::vector<float>& data) {
    GL_CALL(glBindTexture(GL_TEXTURE_1D, texture.id));

    // GL_CALL(glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB32F, data.size() / 3, 0, GL_RGB32F, GL_FLOAT, static_cast<const void*>(data.data()) ));
    GL_CALL(glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB32F, data.size() / 3, 0, GL_RGB, GL_FLOAT, data.data()));

    // if (data.size() != texture.size) {
    //     GL_CALL(glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB32F, data.size(), 0, GL_RGB32F, GL_FLOAT, data.data()));
    // }
    // else {
    //     GL_CALL(glTexSubImage1D(GL_TEXTURE_1D, 0, 0, data.size(), GL_RGB32F, GL_FLOAT, data.data()));
    // }

    // No interpolation in the array
    GL_CALL(glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST));
    GL_CALL(glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST));

    GL_CALL(glBindTexture(GL_TEXTURE_1D, 0));
    texture.size = data.size();
}

LinesGPU::LinesGPU(const std::vector<float>& vertices) {
    // Create VAO and VBO and allocate memory in GPU
    GL_CALL(glGenVertexArrays(1, &VAO));
    GL_CALL(glBindVertexArray(VAO));
    GL_CALL(glGenBuffers(1, &VBO));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, VBO));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW));

    // Set VBO attributes
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, VBO));
    const bool normalized = false;
    GL_CALL(glVertexAttribPointer(0, 3, GL_FLOAT, normalized, 3*sizeof(float), NULL));
    GL_CALL(glEnableVertexAttribArray(0));

    // Disable VBO and VAO
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
    GL_CALL(glBindVertexArray(0));
}

void LinesGPU::drawLines(Material material, const std::vector<float>& vertices) {
    rlEnableShader(material.shader.id);
    // Upload to shader material.colDiffuse
    if (material.shader.locs[SHADER_LOC_COLOR_DIFFUSE] != -1)
    {
        float values[4] = {
            (float)material.maps[MATERIAL_MAP_DIFFUSE].color.r/255.0f,
            (float)material.maps[MATERIAL_MAP_DIFFUSE].color.g/255.0f,
            (float)material.maps[MATERIAL_MAP_DIFFUSE].color.b/255.0f,
            (float)material.maps[MATERIAL_MAP_DIFFUSE].color.a/255.0f
        };

        rlSetUniform(material.shader.locs[SHADER_LOC_COLOR_DIFFUSE], values, SHADER_UNIFORM_VEC4, 1);
    }

    // MATRIX UNIFORMS
    Matrix matModel = MatrixIdentity();
    Matrix matView = rlGetMatrixModelview();
    Matrix matModelView = MatrixIdentity();
    Matrix matProjection = rlGetMatrixProjection();
    Matrix matModelViewProjection = MatrixMultiply(matView, matProjection);

    // Upload view and projection matrices (if locations available)
    if (material.shader.locs[SHADER_LOC_MATRIX_VIEW] != -1) rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_VIEW], matView);
    if (material.shader.locs[SHADER_LOC_MATRIX_PROJECTION] != -1) rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_PROJECTION], matProjection);
    // Model transformation matrix is sent to shader uniform location: SHADER_LOC_MATRIX_MODEL
    if (material.shader.locs[SHADER_LOC_MATRIX_MODEL] != -1) rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_MODEL], matModel);
    // Send combined model-view-projection matrix to shader
    rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_MVP], matModelViewProjection);

    // Update vertices data
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, VBO));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW));

    // Bind mesh VBO data: vertex position (shader-location = 0)
    rlEnableVertexBuffer(VBO);
    rlSetVertexAttribute(material.shader.locs[SHADER_LOC_VERTEX_POSITION], 3, RL_FLOAT, 0, 0, 0);
    rlEnableVertexAttribute(material.shader.locs[SHADER_LOC_VERTEX_POSITION]);

    // Draw the Lines
    GL_CALL(glBindVertexArray(VAO));
    GL_CALL(glDrawArrays(GL_LINES, 0, vertices.size() / 3));

    rlDisableVertexArray();
    rlDisableVertexBuffer();
    rlDisableShader();
    rlSetMatrixModelview(matView);
    rlSetMatrixProjection(matProjection);
}

LinesGPU::~LinesGPU() {
    // Destroy GPU resources
    GL_CALL(glDeleteBuffers(1, &VBO));
    GL_CALL(glDeleteVertexArrays(1, &VAO));
}

inline std::vector<float> static_cast_vector_to_float(const std::vector<Scalar>& vec) {
    std::vector<float> out;
    out.resize(vec.size());
    for (unsigned int i = 0; i < vec.size(); i++)
        out[i] = static_cast<float>(vec[i]);
    return out;
}

SkinnedRodGPU::SkinnedRodGPU(const RenderMesh& mesh, const RodHandle& rod)
{
    Initialize(mesh, rod.bounds.n_rb - 1);
}

SkinnedRodGPU::SkinnedRodGPU(const RenderMesh& mesh, const unsigned int segments) {
    Initialize(mesh, segments);
}

void SkinnedRodGPU::Initialize(const RenderMesh& mesh, const unsigned int segments) {
    this->mesh = mesh;
    // Compute skinning weights
    Eigen::Matrix<Scalar, 2, Eigen::Dynamic> weightsMat;
    Eigen::Matrix<unsigned int, 2, Eigen::Dynamic> boneIDsMat;
    compute_rod_skinning_weights(mesh.vertices, segments, weightsMat, boneIDsMat);

    std::vector<float> weights;
    std::vector<float> boneIDs;
    for (unsigned int i = 0; i < weightsMat.cols(); i++) {
        weights.push_back(static_cast<float>(weightsMat(0,i)));
        weights.push_back(static_cast<float>(weightsMat(1,i)));
        // printf("%f, ", weights[weights.size()-2]);
        // printf("%f\n", weights[weights.size()-1]);

        boneIDs.push_back(static_cast<float>(boneIDsMat(0,i)));
        boneIDs.push_back(static_cast<float>(boneIDsMat(1,i)));
        // printf("%i, ", boneIDs[boneIDs.size()-2]);
        // printf("%i\n", boneIDs[boneIDs.size()-1]);
    }

    // In GPU we want to use floats
    std::vector<float> vertices = static_cast_vector_to_float(mesh.vertices);
    std::vector<float> normals = static_cast_vector_to_float(mesh.normals);
    std::vector<float> tangents = static_cast_vector_to_float(mesh.tangents);
    std::vector<float> texcoords = static_cast_vector_to_float(mesh.texcoords);

    // Create VAO and VBO and allocate memory in GPU
    GL_CALL(glGenVertexArrays(1, &VAO));
    GL_CALL(glBindVertexArray(VAO));

    // Vertices postions
    GL_CALL(glGenBuffers(1, &posVBO));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, posVBO));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW));

    GL_CALL(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), NULL));
    GL_CALL(glEnableVertexAttribArray(0));

    // Vertices normals
    GL_CALL(glGenBuffers(1, &normalVBO));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, normalVBO));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_STATIC_DRAW));

    GL_CALL(glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), NULL));
    GL_CALL(glEnableVertexAttribArray(1));

    // Vertices texture coordinates
    GL_CALL(glGenBuffers(1, &uvVBO));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, uvVBO));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, texcoords.size() * sizeof(float), texcoords.data(), GL_STATIC_DRAW));

    GL_CALL(glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), NULL));
    GL_CALL(glEnableVertexAttribArray(2));

    // Vertices tangents
    GL_CALL(glGenBuffers(1, &tangentVBO));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, tangentVBO));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, tangents.size() * sizeof(float), tangents.data(), GL_STATIC_DRAW));

    GL_CALL(glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), NULL));
    GL_CALL(glEnableVertexAttribArray(3));

    // Vertices bone indices
    GL_CALL(glGenBuffers(1, &idVBO));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, idVBO));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, boneIDs.size() * sizeof(float), boneIDs.data(), GL_STATIC_DRAW));

    GL_CALL(glVertexAttribPointer(4, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), NULL));
    GL_CALL(glEnableVertexAttribArray(4));

    // Vertices bone weights
    GL_CALL(glGenBuffers(1, &weightVBO));
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, weightVBO));
    GL_CALL(glBufferData(GL_ARRAY_BUFFER, weights.size() * sizeof(float), weights.data(), GL_STATIC_DRAW));

    GL_CALL(glVertexAttribPointer(5, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), NULL));
    GL_CALL(glEnableVertexAttribArray(5));

    // Disable VBO and VAO
    GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
    GL_CALL(glBindVertexArray(0));
}

SkinnedRodGPU::~SkinnedRodGPU() {
    // Destroy GPU resources
    GL_CALL(glDeleteBuffers(1, &posVBO));
    GL_CALL(glDeleteBuffers(1, &normalVBO));
    GL_CALL(glDeleteBuffers(1, &uvVBO));
    GL_CALL(glDeleteBuffers(1, &tangentVBO));
    GL_CALL(glDeleteBuffers(1, &idVBO));
    GL_CALL(glDeleteBuffers(1, &weightVBO));
    GL_CALL(glDeleteVertexArrays(1, &VAO));
}

void DrawSkinnedRodGPU(const SkinnedRodGPU& rod, Material material, const std::vector<Matrix>& bone_transforms) {
    rlEnableShader(material.shader.id);
    // Upload to shader material.colDiffuse
    if (material.shader.locs[SHADER_LOC_COLOR_DIFFUSE] != -1)
    {
        float values[4] = {
            (float)material.maps[MATERIAL_MAP_DIFFUSE].color.r/255.0f,
            (float)material.maps[MATERIAL_MAP_DIFFUSE].color.g/255.0f,
            (float)material.maps[MATERIAL_MAP_DIFFUSE].color.b/255.0f,
            (float)material.maps[MATERIAL_MAP_DIFFUSE].color.a/255.0f
        };

        rlSetUniform(material.shader.locs[SHADER_LOC_COLOR_DIFFUSE], values, SHADER_UNIFORM_VEC4, 1);
    }

    // MATRIX UNIFORMS
    Matrix matModel = MatrixIdentity();
    Matrix matView = rlGetMatrixModelview();
    Matrix matProjection = rlGetMatrixProjection();
    Matrix matModelViewProjection = MatrixMultiply(matView, matProjection);

    // Upload view and projection matrices (if locations available)
    if (material.shader.locs[SHADER_LOC_MATRIX_VIEW] != -1) rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_VIEW], matView);
    if (material.shader.locs[SHADER_LOC_MATRIX_PROJECTION] != -1) rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_PROJECTION], matProjection);
    // Model transformation matrix is sent to shader uniform location: SHADER_LOC_MATRIX_MODEL
    if (material.shader.locs[SHADER_LOC_MATRIX_MODEL] != -1) rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_MODEL], matModel);
    if (material.shader.locs[SHADER_LOC_MATRIX_NORMAL] != -1) rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_NORMAL], MatrixTranspose(MatrixInvert(matModel)));
    // Send combined model-view-projection matrix to shader
    rlSetUniformMatrix(material.shader.locs[SHADER_LOC_MATRIX_MVP], matModelViewProjection);

    const unsigned int boneMatricesID = rlGetLocationUniform(material.shader.id, "bonesMatrices");
    std::vector<float> matrices;
    const unsigned int n_transforms = bone_transforms.size();
    for (unsigned int i = 0; i < n_transforms; i++) {
        Matrix mat = bone_transforms[i];
        matrices.insert(matrices.end(), {
            mat.m0, mat.m1, mat.m2, mat.m3,
            mat.m4, mat.m5, mat.m6, mat.m7,
            mat.m8, mat.m9, mat.m10, mat.m11,
            mat.m12, mat.m13, mat.m14, mat.m15
        });
    }
    GL_CALL(glUniformMatrix4fv(boneMatricesID, n_transforms, GL_FALSE, matrices.data()));

    // Draw the Lines
    const unsigned int n_triangles = rod.mesh.vertices.size() / 3;
    GL_CALL(glBindVertexArray(rod.VAO));
    GL_CALL(glDrawArrays(GL_TRIANGLES, 0, n_triangles));

    rlDisableVertexArray();
    rlDisableVertexBuffer();
    rlDisableShader();
    rlSetMatrixModelview(matView);
    rlSetMatrixProjection(matProjection);
}
