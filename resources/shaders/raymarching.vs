// -*- glsl -*-
#version 330

// Input vertex attributes
in vec3 vertexPosition;

void main()
{
    gl_Position = vec4(vertexPosition.xy, 0.0, 1.0);
}
