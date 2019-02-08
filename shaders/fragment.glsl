#version 330

out vec4 colour;
uniform float time;
uniform vec2 resolution;
// Define a sdf of a single sphere
float EPSILON = 0.004;
int MAX_N_STEPS = 200;

// Sphere centered in the origin
// with radius 's'
float sdSphere(vec3 p, float r){
    return length(p) - r;
}
float SDF(vec3 position){
    float t = sdSphere(position - vec3(0.0, 0.0, 10.0), 3.0);
    return t;
}

float RayMarching(vec3 rayOrigin, vec3 rayDirection){
    float t = 0.0;
    for(int i = 0; i < MAX_N_STEPS; i++){
        float dist = SDF(rayOrigin + rayDirection*t);
        if(dist < 0.0001*t){
            return t;
        }
        t += dist;
    }
    return -1.0;
}

vec3 render(vec3 rayOrigin, vec3 rayDirection){
    float t = RayMarching(rayOrigin, rayDirection);
    vec3 col = vec3(1.0 -t*0.075);
    return col;
}
vec2 normalizeScreenCoords(vec2 screenCoord){
    vec2 result = 2.0 * (screenCoord/resolution.xy - 0.5);
    result.x *= resolution.x/resolution.y; // Correct for aspect ratio
    return result;
}
vec3 getCameraRayDir(vec2 uv, vec3 camPos, vec3 camTarget)
{
    // Calculate camera's "orthonormal basis", i.e. its transform matrix components
    vec3 camForward = normalize(camTarget - camPos);
    vec3 camRight = normalize(cross(vec3(0.0, 1.0, 0.0), camForward));
    vec3 camUp = normalize(cross(camForward, camRight));
     
    float fPersp = 2.0;
    vec3 vDir = normalize(uv.x * camRight + uv.y * camUp + camForward * fPersp);
 
    return vDir;
}

void main(){
    // Compute the camera origin
    vec3 camPos = vec3(0, 0, -1);
    vec3 camTarget = vec3(0, 0, 0);
    
    vec2 uv = normalizeScreenCoords(gl_PointCoord);

    // Compute the ray direction
    vec3 rayDir = getCameraRayDir(uv, camPos, camTarget);
    // Compute distance & colour to one surface Raymarching
    vec3 col = render(camPos, rayDir);

    colour = vec4(0.5,0.2,0.2, 1);
}