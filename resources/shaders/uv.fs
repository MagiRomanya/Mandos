// -*- glsl -*-
#version 330

// Input vertex attributes (from vertex shader)
in vec2 fragTexCoord;

// Output fragment color
out vec4 finalColor;


void main()
{
    finalColor = vec4(fragTexCoord, 1.0, 1.0);
}
