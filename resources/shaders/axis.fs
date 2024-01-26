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

#define     LIGHT_DIRECTIONAL       0
#define     LIGHT_POINT             1

struct Light {
    int enabled;
    int type;
    vec3 position;
    vec3 target;
    vec4 color;
};

// Input lighting values
vec4 ambient = vec4(0.2);
vec3 viewPos = vec3(0.0, 0.0, 1.0);

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

    // Texel color fetching from texture sampler
    vec3 lightDot = vec3(0.0);
    vec3 normal = normalize(fragNormal);
    vec3 viewD = normalize(viewPos - fragPosition);
    vec3 specular = vec3(0.0);

    Light light;

    light.enabled = 1;
    light.type = LIGHT_DIRECTIONAL;
    light.position = vec3(0, 0, 4);
    light.target = vec3(0, 0, 0);
    light.color = vec4(1); // white light

    if (light.enabled == 1){
      vec3 light_vec = vec3(0.0);

      if (light.type == LIGHT_DIRECTIONAL)
      {
        light_vec = normalize(light.position - light.target);
      }

      if (light.type == LIGHT_POINT)
      {
        light_vec = normalize(light.position - fragPosition);
      }

      float NdotL = max(dot(normal, light_vec), 0.0);
      lightDot += light.color.rgb * NdotL;

      float specCo = 0.0;
      if (NdotL > 0.0) specCo = pow(max(0.0, dot(viewD, reflect(-(light_vec), normal))), 16.0); // 16 refers to shine
      specular += specCo;
    }

    finalColor = (modelColor + vec4(specular, 1.0))*vec4(lightDot, 1.0);
    finalColor += (ambient)*modelColor;

    // Gamma correction
    finalColor = pow(finalColor, vec4(1.0/2.2));
}
