#version 330

out vec4 color;
in vec2 screenPos;
uniform float Time;
uniform vec2 resolution;
// Camera Parameters
uniform vec3 cameraPosition;
uniform vec3 cameraFront;
uniform vec3 cameraUp;
uniform vec3 cameraRight;
uniform vec3 cameraWorldUp;
uniform mat4 model;
uniform mat4 view;
uniform float screenWidth;
uniform float screenHeight;
uniform vec2 MouseXY;
// Variable for sample textures
uniform sampler3D lowFrequencyTexture;
in vec3 cameraP;
// Define a sdf of a single sphere
float EPSILON = 0.004;
int MAX_N_STEPS = 200;
// Camera Parameters
float CameraSpeed = 5.0;


// Sphere centered in the origin
// with radius 's'
float sdSphere(vec3 p, float r){
    return length(p) - r;
}
float sdPlane( vec3 p){
    return p.y;
}
float SDF(vec3 position){
    float t = sdSphere(position - vec3(0.0, 0.0, -100.0), 2.0);
    //float t = sdPlane(position);
    return t;
}
vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        SDF(vec3(p.x + EPSILON, p.y, p.z)) -  SDF(vec3(p.x - EPSILON, p.y, p.z)),
        SDF(vec3(p.x, p.y + EPSILON, p.z)) -  SDF(vec3(p.x, p.y - EPSILON, p.z)),
        SDF(vec3(p.x, p.y, p.z  + EPSILON)) - SDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

vec4 RayMarching(vec3 rayOrigin, vec3 rayDirection){
    float t = 0.0;
    for(int i = 0; i < MAX_N_STEPS; i++){
        float dist = SDF(rayOrigin + rayDirection*t);
        if(dist < 0.0001*t){
            vec3 normal = estimateNormal(rayOrigin + rayDirection*t);
            return vec4(normal*0.5 + 0.5, t);
            //return t;
        }
        t += dist;
    }
    return vec4( vec3(0.0), -1.0);
}

vec3 render(vec3 rayOrigin, vec3 rayDirection){
    //float t = RayMarching(rayOrigin, rayDirection);
    vec4 res = RayMarching(rayOrigin, rayDirection);
    //vec3 col = vec3(1.0 -t*0.075);
    //float smp = texture(lowFrequencyTexture, vec3(0.2, 0.2, 0.2)).x;
    //vec3 dd = vec3()
    vec4 smp = texture(lowFrequencyTexture, res.xyz);
    vec3 col;// = smp * res.yzw;
    if(res.w < 0.0){
        col = vec3(0.0);
    }else{
        col = normalize(smp.xyz);
    }
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

mat3 setCamera(vec3 ro, vec3 ta, float cr )
{
	vec3 cw = normalize(ta-ro);
	vec3 cp = vec3(sin(cr), cos(cr),0.0);
	vec3 cu = normalize( cross(cw,cp) );
	vec3 cv = normalize( cross(cu,cw) );
    return mat3( cu, cv, cw );
}
mat4 viewMatrix(vec3 front, vec3 right, vec3 up) {
	/*vec3 f = normalize(center - eye);
	vec3 s = normalize(cross(f, up));
	vec3 u = cross(s, f);
	return mat4(
		vec4(s, 0.0),
		vec4(u, 0.0),
		vec4(-f, 0.0),
		vec4(0.0, 0.0, 0.0, 1)
	);*/
    return mat4(
        vec4(right, 0.0),
        vec4(up, 0.0),
        vec4(-front, 0.0),
        vec4(0.0, 0.0, 0.0, 1.0)
    );
}
void main(){
    // Compute the camera origin
    //vec3 camPos = vec3(0, 0, -0.01);
    float aspectratio = screenWidth/screenHeight;
    //camPos += vec3(sin(0.5*Time)*0.5, cos(0.5*Time)*0.1, 0.0);
    //vec3 camPos = model * vec4(0.0, 0.0, -5.0, 1.0);
    vec3 camPos = cameraPosition;// - vec3(0 , 0, -2.0);
    // vec3 camTarget = vec3(0, 0, -1);
    vec2 m = vec2(MouseXY.x/screenWidth, MouseXY.y/screenHeight);
    vec3 r0 = 4.0*normalize(vec3(sin(3.0*m.x), 0.4*m.y, cos(3.0*m.x)));
    float x = aspectratio * (2.0 * (gl_FragCoord.x) / screenWidth - 1.0);
    float y = 2.0 * (gl_FragCoord.y) / screenHeight - 1.0;
    //vec3 v0 = view * vec3(x, y, 0.0);
    vec3 v0 = (view * vec4(x, y, 0.0, 1.0)).xyz;
    //vec3 v0 = (viewMatrix(cameraFront, cameraRight, cameraUp) * vec4(x, y, camPos.z + 0.2, 1.0)).xyz;
    //vec3 rayDir = normalize(v0 - camPos);
    vec3 rayDir = 10.0*cameraFront + cameraRight * (x) + cameraUp * (y); 
    rayDir = normalize(rayDir);
    //vec4 rr = view * model * vec4(rayDir, 0.0);
    mat3 ca = setCamera(r0, v0, 0.0);
    //vec3 rd = ca * normalize(vec3(-vec2(screenWidth, screenHeight)/screenHeight + 2.0*gl_FragCoord.xy/screenHeight, 1.5));

    //vec2 uv = normalizeScreenCoords(gl_FragCoord.xy);

    // Compute the ray direction
    //vec3 rayDir = getCameraRayDir(uv, camPos, camTarget);
    // Compute distance & colour to one surface Raymarching
    vec3 col = render(camPos, rayDir);

    color = vec4(col, 1.0);
}