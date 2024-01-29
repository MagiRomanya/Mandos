// -*- glsl -*-
#version 330

uniform vec4 colDiffuse;

// Output fragment color
out vec4 finalColor;

void main()
{
    finalColor = colDiffuse;
    if (!gl_FrontFacing) {
        finalColor = vec4(vec3(1.0) - finalColor.xyz, 1.0);
    }
}
