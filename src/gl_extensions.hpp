#ifndef GL_EXTENSIONS_H_
#define GL_EXTENSIONS_H_

#include "raylib.h"
#include <vector>

/**
 * Creates a FBO with multisampled color and depth textures.
 */
void InitializeMultisampleFramebuffer();

void BindMultisampleFramebuffer();

/**
 * After all has been rendered onto the Multisampled FBO, all the samples must be combined
 * into a single one to a separate RenderTexture.
 *
 * https://learnopengl.com/Advanced-OpenGL/Anti-Aliasing
 */
void BlitMultisampleFramebuffer(const RenderTexture destination);

void EnableUserDefinedClipping();

void DisableUserDefinedClipping();

struct Texture1D {
    unsigned int id;
    void* data = nullptr;
    int size = 0;
};

Texture1D CreateTexture1D();

void UpdateTexture1D(Texture1D texture, const std::vector<float>& data);

void UseTexture1D(Texture1D& texture);

/**
 * Updates the color and depth attachments resolutions of the FBO.
 * This function also updates the MultiSampled FBO resolution if it is in use.
 */
void UpdateRenderTexture2D(RenderTexture2D& fbo, int width, int height);

struct LinesGPU {
    unsigned int VAO, VBO;

    LinesGPU(const std::vector<float>& vertices);

    void drawLines(Material material, const std::vector<float>& vertices);

    ~LinesGPU();
};

#endif // GL_EXTENSIONS_H_
