// -*- glsl -*-
#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;

uniform vec2 Resolution;
uniform vec3 viewPos;
uniform vec3 viewTarget;
uniform mat4 matView;

// Output fragment color
out vec4 finalColor;

float sdCapsule( vec3 p, vec3 a, vec3 b, float r ) {
  vec3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// Scene distance
float map(vec3 p) {
    // return sdBox(p, vec3(1.0,1.0,1.0));
    return length(p) - 1.0; // distance to a sphere of radius 1
    // float result = sdCapsule(p, vec3(-1,0,0), vec3(-0.8,0.50,0.2), 0.2);

    // return result;
}
vec3 calcNormal( in vec3 pos )
{
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*map( pos + e.xyy ) +
                      e.yyx*map( pos + e.yyx ) +
                      e.yxy*map( pos + e.yxy ) +
                      e.xxx*map( pos + e.xxx ) );
}


mat3 setCamera( in vec3 ro, in vec3 ta, float cr )
{
    vec3 cw = normalize(ta-ro);
    vec3 cp = vec3(sin(cr), cos(cr), 0.0);
    vec3 cu = normalize( cross(cw,cp) );
    vec3 cv = normalize( cross(cu,cw) );
    return mat3( cu, cv, cw );
}

void main() {
    vec2 uv = (gl_FragCoord.xy * 2.0 - Resolution.xy) / Resolution.y;

    // Initialization
    vec3 ro = viewPos;
    vec3 ta = viewTarget;
    mat3 ca = setCamera( ro, ta, 0.0 );
    // ray direction
    float fov = radians(45.0);  // convert degrees to radians
    vec2 uv_fov = vec2(uv.x * tan(fov * 0.5), uv.y * tan(fov * 0.5));
    // vec3 rd = ca * normalize( vec3(uv, 2.0) );
    vec3 rd = ca * normalize( vec3(uv_fov, 1.0) );
    vec3 col = vec3(0);               // final pixel color

    float t = 0.0; // total distance travelled

    // Raymarching
    for (int i = 0; i < 80; i++) {
        vec3 p = ro + rd * t;     // position along the ray

        float d = map(p);         // current distance to the scene

        t += d;                   // "march" the ray

        // col = vec3(i) / 80.0;
        if (d < 0.001) break;      // early stop if close enough
        if (t > 100.0) break;      // early stop if too far
    }
    float near=  0.01;
    float far = 1000.0;
    gl_FragDepth = (1/t  -1/near) / (1/far - 1/near);

    if (t > 100.0) {
        finalColor = vec4(0.0, 0.0, 0.0, 0.0);
        return;
    }

    // Coloring
    col = vec3(t * 0.2);           // color based on distance

    finalColor = vec4(col, 1);
}
