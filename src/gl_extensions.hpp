#ifndef GL_EXTENSIONS_H_
#define GL_EXTENSIONS_H_

#include "raylib.h"
#include <vector>

struct LinesGPU {
    unsigned int VAO, VBO;

    LinesGPU(const std::vector<float>& vertices);

    void drawLines(Material material, const std::vector<float>& vertices);

    ~LinesGPU();
};

#endif // GL_EXTENSIONS_H_
