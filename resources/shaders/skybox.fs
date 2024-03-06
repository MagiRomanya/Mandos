// -*- glsl -*-
#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;

// Input uniform values
uniform samplerCube environmentMap;

// Output fragment color
out vec4 finalColor;

void main()
{
    // Fetch color from texture map
    vec3 color = texture(environmentMap, fragPosition).rgb;

    finalColor = vec4(color, 1.0);
}
