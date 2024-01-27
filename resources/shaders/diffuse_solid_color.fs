// -*- glsl -*-
#version 330

uniform vec4 colDiffuse;

// Output fragment color
out vec4 finalColor;

void main()
{
    finalColor = colDiffuse;
}
