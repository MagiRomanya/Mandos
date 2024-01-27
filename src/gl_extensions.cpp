#include <cassert>
#include <iostream>
#include <vector>

#include "external/glad.h"
#include "gl_extensions.hpp"
#include "raymath.h"
#include <rlgl.h>

// OPENGL ERROR HANDLING
// ----------------------------------------------------------------
static void GLClearError() {
    while (glGetError() != GL_NO_ERROR);
}

static bool GLLogCall(char* call_name, char* file, int line) {
    while (GLenum error = glGetError()) {
        std::cerr << "[OpenGL ERROR] (" << error << ") in file "
            << file << ", and line " << line <<  std::endl;
        return false;
    }
    return true;
}

#define GL_CALL(x) \
    GLClearError(); \
    x; \
    assert(GLLogCall(#x, __FILE__, __LINE__))
// ----------------------------------------------------------------


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