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
uniform sampler2D WeatherTexture;
uniform sampler1D GradientStratusTexture;
uniform sampler1D GradientCumulusTexture;
uniform sampler1D GradientCumulonimbusTexture;

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

// Remap from book

float Remap(float original_value, float original_min, float original_max, float new_min, float new_max){
    return new_min + (((original_value - original_min)/(original_max - original_min)) * (new_max - new_min));
}


/* @brief
 * Density height function based on wheather map
 */
float GetDensityHeightGradientForPoint(vec2 point, float height){
    // Sample the texture pn point
    vec4 weatherdata = texture(WeatherTexture, point);
    // Cloud coverage on Sky this is red channel
    float cloudcoverage = weatherdata.x;
    // precipitation probability
    float precipitation = weatherdata.y;
    // cloud type {0: Stratus, 1: stratocummulus, 2: cumulus}
    int cloudtype = 0;
    if(weatherdata.z < 0.55 && weatherdata.z > 0.45)
        cloudtype = 1;
    else if(weatherdata.z > 0.9)
        cloudtype = 2;

    // TODO: for moment just cloud coverage is been using.
    return cloudcoverage;

}

/*
 * @brief
 * Get the texture weather map intesection
 * by mapping the texture over all the dome
 */
vec2 GetPointToSampleInDome(vec3 positionInDome, vec3 eyePosition, vec3 domeCenter, float radiusDome){
    float r = length(eyePosition - domeCenter);
    float dist = sqrt(radiusDome*radiusDome - r * r);
    vec3 posXmin = eyePosition + vec3(-1.0, 0.0, 0.0) * dist;
    vec3 posZmin = eyePosition + vec3(0.0, 0.0, -1.0) * dist;
    
    posXmin.y = 0.0;
    posZmin.y = 0.0;
    positionInDome.y = 0.0;


    float posx = length(positionInDome - posXmin)/(2 * dist);
    float posy = length(positionInDome - posZmin)/(2 * dist);
    vec2 pointSample = vec2(posx, posy);
    return pointSample;
}

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
vec3 GetPositionInAtmosphere(vec3 pos, vec3 earthCenter, float thickness){
    vec3 posit = (pos - vec3(earthCenter.x, ATMOSPHERE_INNER_RADIUS - EARTH_RADIUS, earthCenter.z))/thickness;

    return posit;
}
/* 
 * @brief:
 * Sample Low Frequency Texture - Provides Brownian noise -> FBM
 */

float sampleCloudDensity(vec3 p, vec2 weather_point, float relativeHeight){
    // Read the low frequency Perlin-Worley & Worley noise
    vec4 low_frequency_noises = texture(lowFrequencyTexture, p);
    float low_freq_FBM = low_frequency_noises.y * 0.625 + 
                        low_frequency_noises.z * 0.250 +
                        low_frequency_noises.w * 0.125;
    float base_cloud = Remap(low_frequency_noises.x, -(1.0 - low_freq_FBM), 1.0, 0.0, 1.0);

    // TODO: Base on whether texture
    
    float density_height_gradient = GetDensityHeightGradientForPoint(weather_point, relativeHeight);

    return base_cloud * density_height_gradient; 
}
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
        vec3 pointSampled = GetPositionInAtmosphere(pos, earthCenter, atmosphereThickness);//vec3(0.6, 0.0, 0.6);

        vec2 weatherPoint = GetPointToSampleInDome(pos, rayOrigin, earthCenter, ATMOSPHERE_OUTER_RADIUS);
        float baseCloud = sampleCloudDensity(pointSampled, weatherPoint, relativeHeight);
        //float baseCloud = sampleLowFrequencyTexture(pointSampled);
        //float baseCloud = 0.5;
        if(baseCloud > 0.0){
            acummDensity += baseCloud *0.0001;
            //acummDensity += baseCloud * 0.1;
            colorPixel += vec3(baseCloud *0.006);
            //break;
        }
        if(acummDensity >= 1.0){
            acummDensity = 1.0;
            break;
        }
        if(baseCloud < 0){
            colorPixel = vec3(0.0);
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
        vec3 colorNearHorizon = vec3(0.54, 0.23, 0.046) * 0.4;
        color =  vec4(colorNearHorizon, 1.0);
        return;
    }
    /*else if(seeingUp < 0.06){
        vec3 colorHorizon = vec3(0.3, 0.32, 0.51) * 0.4;
        color =  vec4(colorHorizon, 1.0);
        return;
    }*/
    //vec4 weatherdata = texture(WeatherTexture, vec2(x,y));
    vec4 gradient = texture(GradientCumulonimbusTexture, y);
    vec3 col = vec3(gradient.x);
    //float typecloud = weatherdata.z;
    //float red = 0.0;
    //float green = 0.0;
    //float blue = 0.0;
    //if(weatherdata.z < 0.55 && weatherdata.z > 0.45)
    //    green = 1.0;
    //else if(weatherdata.z > 0.9)
    //    blue = 1.0;
    //else if(weatherdata.z < 0.1)
    //    red = 1.0;
//
    //vec3 col = vec3(red, green, blue);
    //vec3 col = RayMarching(rayOrigin, rayDirection, innerIntersection, earthCenter, initialLength, finalLength);
    //col = vec3(1.0 - col.x, 1.0 - col.y, 1.0 - col.z);
    //float modcolor = length(col);
    //if(modcolor < 0.4)
    //    col = vec3(0.5273, 0.808, 0.922);
    color = vec4(col, 1.0);
}
