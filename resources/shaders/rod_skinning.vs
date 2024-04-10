// -*- glsl -*-
#version 430 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 tex;
layout(location = 3) in vec3 tangent;
layout(location = 4) in vec2 boneIds;
layout(location = 5) in vec2 weights;

uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;

const int MAX_BONES = 200;
const int MAX_BONE_INFLUENCE = 2;
uniform mat4 bonesMatrices[MAX_BONES];

out vec3 fragPosition;
out vec2 fragTexCoord;
out vec3 fragNormal;
out vec3 fragTangent;
out mat3 TBN;

void main()
{
    vec4 totalPosition = vec4(0.0f);
    vec3 totalNormal = vec3(0.0);
    for(int i = 0 ; i < MAX_BONE_INFLUENCE ; i++)
    {
        if(boneIds[i] == -1)
            continue;
        if(boneIds[i] >=MAX_BONES)
        {
            totalPosition = vec4(pos,1.0f);
            break;
        }
        int id = int(boneIds[i]);
        vec4 localPosition = bonesMatrices[id] * vec4(pos,1.0f);
        totalPosition += localPosition * weights[i];
        vec3 localNormal = mat3(bonesMatrices[id]) * normal;
        totalNormal += localNormal * weights[i];
    }

    vec4 projPos = mvp * totalPosition;
    // mat4 scale;
    // scale[0] = vec4(1.0, 0.0, 0.0, 0.0);
    // scale[1] = vec4(0.0, 1.0, 0.0, 0.0);
    // scale[2] = vec4(0.0, 0.0, 1.0, 0.0);
    // scale[3] = vec4(0.0, 0.0, 0.0, 1.0);
    // vec4 projPos = mvp * vec4(pos, 1.0);
    gl_Position =  projPos;


    fragPosition = projPos.xyz;
    fragTexCoord = tex;
    fragTexCoord = boneIds;
    if (int(boneIds[0]) == 0) {
        fragTexCoord = weights;
    }
    else {
        fragTexCoord = weights.yx;
    }
    // fragNormal = normalize(vec3(matNormal*vec4(normal, 1.0)));
    fragNormal = normalize(totalNormal);
    fragTangent = normalize(vec3(matModel*vec4(tangent, 1.0)));
    vec3 T = fragTangent;
    vec3 N = fragNormal;
    vec3 B = normalize(cross(N,T));
    TBN = mat3(T, B, N);
}
