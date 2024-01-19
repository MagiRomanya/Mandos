// -*- glsl -*-
#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec3 fragNormal;
in vec3 modelPosition;

// Input uniform values
uniform sampler2D texture0;
uniform vec4 colDiffuse;

// Output fragment color
out vec4 finalColor;

void main()
{
    vec4 modelColor;
    if (length(modelPosition) > 0.1) {
        float maxCoord = max(max(modelPosition.x, modelPosition.y), modelPosition.z);
        if (maxCoord == modelPosition.x) {
            vec4 red = vec4(1.0, 0.1, 0.1, 1.0);
            modelColor = colDiffuse * red;
        }
        else if (maxCoord == modelPosition.y) {
            vec4 green = vec4(0.1, 1.0, 0.1, 1.0);
            modelColor = colDiffuse * green;
        }
        else if (maxCoord == modelPosition.z) {
            vec4 blue = vec4(0.1, 0.1, 1.0, 1.0);
            modelColor = colDiffuse * blue;
        }
    }
    else {
        modelColor = colDiffuse;
    }
    finalColor = modelColor;
}
