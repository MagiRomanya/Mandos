// -*- glsl -*-
#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;
//in vec4 fragColor;
in vec3 fragNormal;

// Input uniform values
uniform sampler2D texture0;
uniform vec4 colDiffuse;

// Output fragment color
out vec4 finalColor;

// NOTE: Add here your custom variables

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
uniform vec3 viewPos;

void main()
{
    // Texel color fetching from texture sampler
    vec4 texelColor = texture(texture0, fragTexCoord);
    vec3 lightDot = vec3(0.0);
    vec3 normal = normalize(fragNormal);
    vec3 viewD = normalize(viewPos - fragPosition);
    vec3 specular = vec3(0.0);

    // NOTE: Implement here your fragment shader code

    Light light;

    light.enabled = 1;
    light.type = LIGHT_POINT;
    light.position = vec3(0, 8, 4);
    light.target = vec3(0, 0, 0);
    light.color = vec4(1.0); // white light

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

    finalColor = (texelColor*((colDiffuse + vec4(specular, 1.0))*vec4(lightDot, 1.0)));
    finalColor += texelColor*(ambient/10.0)*colDiffuse;

    // Gamma correction
    finalColor = pow(finalColor, vec4(1.0/2.2));
    // finalColor = vec4(normal, 1.0);
}
