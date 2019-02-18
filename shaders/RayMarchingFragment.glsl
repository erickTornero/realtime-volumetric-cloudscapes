#version 330

out vec4 color;

uniform vec3 cameraPosition;
uniform vec3 cameraFront;
uniform vec3 cameraUp;
uniform vec3 cameraRight;

uniform float screenWidth;
uniform float screenHeight;

// Variable for sample textures
uniform sampler3D lowFrequencyTexture;

struct Ray{
    vec3 origin;
    vec3 direction;
};
vec3 GetRayDirection(vec3 front, vec3 right, vec3 up, float x, float y){
    vec3 ray = 10 * front + right * x + up * y;
    return normalize(ray);
}
vec3 GetIntersectionRay2Sphere(vec3 rayOrigin, vec3 rayDirection, vec3 sphereCenter, float radius){
    
}
float getHeightInAtmosphere(vec3 pointInAtm, vec3 earthCenter, vec3 intersectionRay2InnerAtm, vec3 rayDirection, vec3 eyePos, float atmThick){
    float distanceCamera2Point = length(pointInAtm - eyePos);
    float distanceCamera2InnerAtm = length(intersectionRay2InnerAtm - eye);
    // Asumming module of Ray direction
    vec3 point2EarthCenter = normalize(pointInAtm - earthCenter);

    float cosThet = dot(rayDirection, point2EarthCenter);
    
    float posInAtm = abs(cosThet * (distanceCamera2Point - distanceCamera2InnerAtm));

    return posInAtm/atmThick;
}

vec4 RayMarching(vec3 rayOrigin, vec3 rayDirection, vec3 innerIntersection, vec3 earthCenter, float start_atm, float end_atm){
    const float atmosphereThickness = end_atm - start_atm;
    //TODO: Check if maxNSteps is Ok!
    const float maxNSteps = 200.0;
    const float incrementAtm = atmosphereThickness/maxNSteps;

    vec3 colorPixel = vec3(0.0);
    // Start Ray marching
    for(float t = start_atm; t < end_atm; t += incrementAtm){
        vec3 colorSample = vec3(0.0);

        // Compute position
        vec3 pos = rayOrigin + t * rayDirection;

        float relativeHeight = getHeightInAtmosphere(pos, earthCenter, xx, rayDirection, rayOrigin, atmosphereThickness);



    }
}