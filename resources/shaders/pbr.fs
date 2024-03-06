// -*- glsl -*-
#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;
in vec3 fragTangent;
in mat3 TBN;

// Input uniform values
uniform sampler2D texture0; // diffuse lighting
uniform sampler2D texture2; // normal map
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
uniform vec3 viewPos;
const float PI = 3.1415926535;

vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;

    float num   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float num   = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return num / denom;
}
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}

vec3 normalMapping() {
    vec3 Normal = normalize(fragNormal);
    Normal = texture(texture2, fragTexCoord).rgb;
    Normal = Normal * 2.0 - 1.0;
    Normal = normalize(TBN * Normal);
    return Normal;
}

void main()
{
    // Normal mapping
    vec3 Normal = normalMapping();

    // Light definition
    Light light;
    light.enabled = 1;
    light.type = LIGHT_DIRECTIONAL;
    light.position = vec3(0.0, 8.0, 5.0);
    light.target = vec3(0.0, 0.0, 0.0);
    light.color = vec4(1.0); // white light

    /// PBR
    //Material definition
    vec3 albedo = texture(texture0, fragTexCoord).rgb * 4;
    // vec3 albedo = colDiffuse.xyz * 4;
    const float metallic = 0.0;
    const float roughness = 0.3;
    const float ambientOclusion = 1.0;

    // Necessary vectors
    vec3 Lo = vec3(0.0);
    vec3 View = normalize(viewPos - fragPosition);
    vec3 Light = normalize(light.position - fragPosition);
    vec3 Half = normalize(View + Light);
    vec3 radiance = light.color.rgb;

    // Computing the specular component
    // Fresnel ( ratio of surface reflection )
    vec3 F0 = vec3(0.04);
    F0      = mix(F0, albedo, metallic);
    vec3 F  = fresnelSchlick(max(dot(Half, View), 0.0), F0);
    // Normal distribution function
    float NDF = DistributionGGX(Normal, Half, roughness);
    // Geometry term ( microfacets self oclusion )
    float G   = GeometrySmith(Normal, View, Light, roughness);

    // Computing the full BRDF
    vec3 numerator    = NDF * G * F;
    float denominator = 4.0 * max(dot(Normal, View), 0.0) * max(dot(Normal, Light), 0.0)  + 0.0001;
    vec3 specular     = numerator / denominator;

    // Getting the diffuse component from the specular and conservation of energy
    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS;
    kD *= 1.0 - metallic;

    // Assembling all into the render equiation
    float NdotL = max(dot(Normal, Light), 0.0);
    Lo += (kD * albedo / PI + specular) * radiance * NdotL;

    // Ambient lighting
    vec3 ambient = vec3(0.05) * albedo * ambientOclusion;
    vec3 color   = ambient + Lo;
    // vec3 color   = Lo;

    // HDR and Tone Mapping
    // color = color / (color + vec3(1.0));
    // color = pow(color, vec3(1.0/2.2));

    finalColor = vec4(color, 1.0);
}
