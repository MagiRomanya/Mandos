// -*- glsl -*-
#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;
in vec4 vertexColor;
in vec3 vertexTangent;

// Input uniform values
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
uniform vec4 slicePlane;

// Output vertex attributes (to fragment shader)
out vec3 fragPosition;
out vec2 fragTexCoord;
out vec4 fragColor;
out vec3 fragNormal;
out vec3 fragTangent;
out mat3 TBN;

out float gl_ClipDistance[1];

void main()
{
    // Send vertex attributes to fragment shader
    fragPosition = vec3(matModel*vec4(vertexPosition, 1.0));
    fragTexCoord = vertexTexCoord;
    fragColor = vertexColor;
    fragNormal = normalize(vec3(matNormal*vec4(vertexNormal, 1.0)));
    fragTangent = normalize(vec3(matModel*vec4(vertexTangent, 1.0)));
    vec3 T = fragTangent;
    vec3 N = fragNormal;
    vec3 B = normalize(cross(N,T));
    TBN = mat3(T, B, N);

    // Calculate final vertex position
    gl_Position = mvp*vec4(vertexPosition, 1.0);

    // Plane culling
    gl_ClipDistance[0] = -dot(slicePlane, vec4(fragPosition, 1.0));
}
