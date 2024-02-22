#include <cassert>
#include <iostream>
#include <vector>

#include "external/glad.h"
#include "gl_extensions.hpp"
#include "raymath.h"
#include "utility_functions.hpp"
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

void UpdateRenderTexture2D(RenderTexture2D& fbo, int width, int height) {
    if (width <= 0 or height <= 0) return;
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

    rlSetFramebufferWidth(width);
    rlSetFramebufferHeight(height);
    fbo.texture.width = width;
    fbo.texture.height = height;
    fbo.depth.width = width;
    fbo.depth.height = height;
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
