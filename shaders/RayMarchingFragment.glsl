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

// ** Initialize Global Variables
// Earth radius in meters
const float EARTH_RADIUS = 6378000.0;
const float ATMOSPHERE_INNER_RADIUS = EARTH_RADIUS + 10000.0;
const float ATMOSPHERE_THICKNESS = 10000.0;
const float ATMOSPHERE_OUTER_RADIUS = ATMOSPHERE_INNER_RADIUS + ATMOSPHERE_THICKNESS;

/*
 *  @brief
 *  remaping
 */
float remapClamped(float value, float old_min, float old_max, float new_min, float new_max){
    float v = new_min +((value - old_min)/ (old_max - old_min))*(new_max - new_min);
    return clamp(v, new_min, new_max);
}
float remapClampPrevPost(float value, float old_min, float old_max, float new_min, float new_max){
    value = clamp(value, old_min, old_max);
    return remapClamped(value, old_min, old_max, new_min, new_max);
}

vec3 GetRayDirection(vec3 front, vec3 right, vec3 up, float x, float y){
    vec3 ray = 10 * front + right * x + up * y;
    return normalize(ray);
}
// Get the intersection with a ray launch from 'rayOrigin', whit a sphere with 'radius'
vec3 GetIntersectionRay2Sphere(vec3 rayOrigin, vec3 rayDirection, vec3 sphereCenter, float radius){
    // Solve the cuadratic equation
    float a = 1.0;
    float b = 2 * dot(rayOrigin - sphereCenter, rayDirection);
    float c = dot(rayOrigin - sphereCenter, rayOrigin - sphereCenter) - radius * radius;
    // Solve second grade equation
    float discriminant = b * b - 4 * a * c;
    float t;
    if(discriminant < 0.0){
        // ** 
        return vec3(-1.0);
    }
    else{
        t = (-b - sqrt(discriminant))/(2.0 * a);
        if(t < 0.0)
            t = (-b + sqrt(discriminant))/(2.0 * a);
    }
    return rayOrigin + t * rayDirection;
}
/*
 *  Get the Relative Height [0, 1] in the considered atmosphere
 */
float GetHeightInAtmosphere(vec3 pointInAtm, vec3 earthCenter, vec3 intersectionRay2InnerAtm, vec3 rayDirection, vec3 eyePos, float atmThick){
    float distanceCamera2Point = length(pointInAtm - eyePos);
    float distanceCamera2InnerAtm = length(intersectionRay2InnerAtm - eyePos);
    // Asumming module of Ray direction
    vec3 point2EarthCenter = normalize(pointInAtm - earthCenter);

    float cosThet = dot(rayDirection, point2EarthCenter);
    
    float posInAtm = abs(cosThet * (distanceCamera2Point - distanceCamera2InnerAtm));

    return posInAtm/atmThick;
}
/*
 * ** Relative position in atmosphere
 */
vec3 GetPositionInAtmosphere(vec3 pos, vec3 earthCenter){
    vec3 posit = (pos - vec3(earthCenter.x, ATMOSPHERE_INNER_RADIUS - EARTH_RADIUS, earthCenter.x))/10000.0;

    return posit;
}
/* 
 * @brief:
 * Sample Low Frequency Texture - Provides Brownian noise -> FBM
 */
float sampleLowFrequencyTexture(vec3 pointSample){
    // 3D texture sampled 4 channels RGBA
    vec4 sampledTexture = texture(lowFrequencyTexture, pointSample);
    float valueFBM = sampledTexture.y * 0.625 + sampledTexture.z * 0.25 + sampledTexture.w * 0.125;

    valueFBM = clamp(valueFBM, 0.0, 1.0);
    float coverage = 0.75;
    float baseCloud = remapClamped(sampledTexture.x, (valueFBM - 0.9), 1.0, 0.0, 1.0);
    float baseCloudCoverage = remapClampPrevPost(baseCloud, coverage, 1.0, 0.0, 1.0);

    baseCloudCoverage *= coverage;

    return baseCloudCoverage;
}

vec3 RayMarching(vec3 rayOrigin, vec3 rayDirection, vec3 innerIntersection, vec3 earthCenter, float start_atm, float end_atm){
    float atmosphereThickness = end_atm - start_atm;
    //TODO: Check if maxNSteps is Ok!
    float maxNSteps = 200.0;
    float incrementAtm = atmosphereThickness/maxNSteps;

    float acummDensity = 0.0;
    vec3 colorPixel = vec3(0.0);
    // Start Ray marching
    for(float t = start_atm; t < end_atm; t += incrementAtm){
        vec3 colorSample = vec3(0.0);

        // Compute position to sample in atmosphere
        vec3 pos = rayOrigin + t * rayDirection;

        float relativeHeight = GetHeightInAtmosphere(pos, earthCenter, innerIntersection, rayDirection, rayOrigin, atmosphereThickness);

        // Warn: Hardcoded variable
        vec3 pointSampled = GetPositionInAtmosphere(pos, earthCenter);//vec3(0.6, 0.0, 0.6);

        float baseCloud = sampleLowFrequencyTexture(pointSampled);
        if(baseCloud > 0.0){
            acummDensity += baseCloud *0.5;
            colorPixel += vec3(acummDensity );
        }
        if(acummDensity >= 1.0){
            acummDensity = 1.0;
            break;
        }

    }
    return colorPixel;
}

void main(){
    float aspecRatio = screenWidth / screenHeight;
    vec3 rayOrigin = cameraPosition;
    float x = aspecRatio * (2.0 * gl_FragCoord.x/screenWidth - 1.0);
    float y = 2.0 * gl_FragCoord.y / screenHeight - 1.0;

    // ** Get normalized ray direction
    vec3 rayDirection = GetRayDirection(cameraFront, cameraRight, cameraUp, x, y);

    // ** Start Variables to Ray Marching
    vec3 earthCenter = vec3(rayOrigin.x, rayOrigin.y - EARTH_RADIUS, rayOrigin.z);
    vec3 innerIntersection = GetIntersectionRay2Sphere(rayOrigin, rayDirection, earthCenter, ATMOSPHERE_INNER_RADIUS);
    vec3 outerIntersection = GetIntersectionRay2Sphere(rayOrigin, rayDirection, earthCenter, ATMOSPHERE_OUTER_RADIUS);

    float initialLength = length(innerIntersection - rayOrigin);
    float finalLength = length(outerIntersection - rayOrigin);

    float seeingUp = dot(vec3(0.0, 1.0, 0.0), rayDirection);
    if(seeingUp < 0.0){
        vec3 colorNearHorizon = vec3(0.0, 0.16, 0.51) * 0.4;
        color =  vec4(colorNearHorizon, 1.0);
        return;
    }
    else if(seeingUp < 0.06){
        vec3 colorHorizon = vec3(0.3, 0.32, 0.51) * 0.4;
        color =  vec4(colorHorizon, 1.0);
        return;
    }
    vec3 col = RayMarching(rayOrigin, rayDirection, innerIntersection, earthCenter, initialLength, finalLength);

    color = vec4(col, 1.0);
}
