#version 330

out vec4 color;

uniform vec3 cameraPosition;
uniform vec3 cameraFront;
uniform vec3 cameraUp;
uniform vec3 cameraRight;

uniform float screenWidth;
uniform float screenHeight;
uniform float Time;

uniform sampler3D lowFrequencyTexture;
uniform sampler3D highFrequencyTexture;
uniform sampler2D WeatherTexture;

uniform sampler2D CurlNoiseTexture;

// Define Global variables used by all the program
const float MPI = 3.14159265;
const float EARTH_RADIUS = 6378000.0;
const float THICK_ATMOSPHERE = 3500.0;
const float ATMOSPHERE_INNER_RADIUS = EARTH_RADIUS + 1500.0;
const float ATMOSPHERE_OUTER_RADIUS = ATMOSPHERE_INNER_RADIUS + THICK_ATMOSPHERE;

// The Sun Location
const vec3 SunLocation = vec3(0.0, ATMOSPHERE_OUTER_RADIUS * 0.9, -ATMOSPHERE_OUTER_RADIUS * 0.9);

// ** Definition of samplers

// 3D Texture - Return the RGBA sample point Perlin - Worley nose
vec4 SampleLowFrequencyTexture(vec3 point){
    return texture(lowFrequencyTexture, point);
}
// Sample High Frequency Texture 3D RGB 
vec3 SampleHighFrequencyTexture(vec3 point){
    return texture(highFrequencyTexture, point).xyz;
}
// Sample Weather Texture 2D RGB 
vec3 SampleWeatherTexture(vec2 point){
    return texture(WeatherTexture, point).xyz;
}
// Sample Curl noise texture 2D RGB
vec3 SampleCurlNoiseTexture(vec2 point){
    return texture(CurlNoiseTexture, point).xyz;
}

// ** Get ray direction definition
vec3 GetRayDirection(vec3 front, vec3 right, vec3 up, float x, float y){
    vec3 ray = 5 * front + right * x + up * y;
    return normalize(ray);
}
// ** Remap function defined by the author of paper
float Remap(float original_value, float original_min, float original_max, float new_min, float new_max){
    return new_min + (((original_value - original_min)/(original_max - original_min)) * (new_max - new_min));
}

// GetHeighfraction

float GetHeightFractionForPoint(vec3 inPosition, vec2 inCloudMinMax){
    float height_fraction = (inPosition.z - inCloudMinMax.x)/(inCloudMinMax.y - inCloudMinMax.x);
    return clamp(height_fraction, 0.0, 1.0);
}

// Keep in box [0, 1]
vec3 KeepInBox(vec3 point){
    vec3 newpoint = point;
    if(point.x < 0)
        newpoint.x = 1.0 - point.x;
    if(point.y < 0)
        newpoint.y = 1.0 - point.y;
    if(point.z < 0)
        newpoint.z = 1.0 - point.z;
    
    return newpoint;
}

// Get Relative cloud in thick atmosphere
float GetRelativeHeightInAtmosphere(vec3 position, vec3 earthCenter){
    float distCenter2Point = length(position - earthCenter);
    float relativeDistance = distCenter2Point - ATMOSPHERE_INNER_RADIUS;

    return relativeDistance/THICK_ATMOSPHERE;
}

// temporal computing of height
float GetHeightInAtmosphere(vec3 pointInAtm, vec3 earthCenter, vec3 intersectionRay2InnerAtm, vec3 rayDirection, vec3 eyePos, float atmThick){
    float distanceCamera2Point = length(pointInAtm - eyePos);
    float distanceCamera2InnerAtm = length(intersectionRay2InnerAtm - eyePos);
    // Asumming module of Ray direction
    vec3 point2EarthCenter = normalize(pointInAtm - earthCenter);

    float cosThet = dot(rayDirection, point2EarthCenter);
    
    float posInAtm = abs(cosThet * (distanceCamera2Point - distanceCamera2InnerAtm));

    return clamp(posInAtm/atmThick, 0.0, 1.0);
}

/*
 * @brief: GetGradientHeightFactor
 * Defines the gradient height functions, up & down type sin.
 * The function parameters is defined by the cloudtype
 * cloudtype:  Type of cloud: 0-> Stratus, 1-> Cumulus, 2-> Cumulonimbus
 * Height:     Relative height [0-1] respect to the atmosphere thickness 
 */

float GetGradientHeightFactor(float height, int cloudtype){
    float timewidthup, starttimeup, starttimedown;
    // Stratus clouds
    if(cloudtype == 0){
        timewidthup     = 0.08;
        starttimeup     = 0.08;
        starttimedown   = 0.2;
    } // Cumulus clouds
    else if(cloudtype == 1){
        timewidthup     = 0.14;
        starttimeup     = 0.1;
        starttimedown   = 0.5;
        }
    // Cumulunimbusclouds
    else if(cloudtype == 2){
        timewidthup     = 0.2;
        starttimeup     = 0.10;
        starttimedown   = 0.7;
    }

    float factor = 2.0 * MPI/(2.0 * timewidthup);

    float density_gradient = 0.0;
    // Gradient functions dependent on the cloudtype
    if(height < starttimeup)
        density_gradient = 0.0;
    else if(height < starttimeup + timewidthup)
        density_gradient = 0.5 * sin(factor * height - MPI/2.0 - factor * starttimeup) + 0.5;
    else if(height < starttimedown)
        density_gradient = 1.0;
    else if(height < starttimedown + timewidthup)
        density_gradient = 0.5 * sin(factor * height - MPI/2.0 - factor * (starttimedown + timewidthup)) + 0.5;
    else
        density_gradient = 0.0;

    return density_gradient;
}

// TODO: Remove argument point if it's not necessary:

/*
 * @brief: GetDensityHeightGradientForPoint
 * Return the gradient in accordance to the cloudtype and the Relative Height.
 * vec3  weather_data:   data vector for the sampling point z: cloudtype
 * float relativeHeight: Relative height respect to the atmosphere thickness
 */
float GetDensityHeightGradientForPoint(vec3 point, vec3 weather_data, float relativeHeight){
    float cloudt = weather_data.z;
    // Cloud type: {0: Stratus, 1: Cumulus, 2: Cumulonimbus}
    int cloudtype = 1;
    if(cloudt < 0.1)
        cloudtype = 0;
    else if(cloudt > 0.9)
        cloudtype = 2;
    
    // Relative Height from [0 - 1]
    //float relh = point.y;
    //float relativeHeight = GetRelativeHeightInAtmosphere(point, earthCenter);

    // Get gradient function defined for three type of clouds
    return GetGradientHeightFactor(relativeHeight, cloudtype);
}
// Define the intersection between ray & dome

vec3 GetIntersectionSphereRay(vec3 rayOrigin, vec3 rayDirection, vec3 sphereCenter, float sphereRadius){
    // Solve an second grade equation
    // Distance from spherecenter to ray origin
    float er = length(rayOrigin - sphereCenter);
    // Vector of direction normalized
    vec3 uer = normalize(rayOrigin - sphereCenter);
    // if Intersection = rayOrigin + rayDirection * t
    // The second order equation is:
    // er^2 + t^2 + er*t*dot(uer, rayDirection) = sphereRadius^2
    float A = 1.0;
    float B = 2 * er * dot(uer, rayDirection);
    float C = er * er - sphereRadius * sphereRadius;
    // Define the discriminant
    float discriminant = B * B - 4 * A * C;
    // Solve for t the equation
    float t;
    if(discriminant < 0.0)
        return vec3(0.0);
    else{
        t = (-B - sqrt(discriminant))/(2.0 * A);
        if(t < 0.0)
            t = (-B + sqrt(discriminant))/(2.0 * A);
    }

    return rayOrigin + t * rayDirection;
}

// Get point sampling 3D texture

vec3 GetPositionInAtmosphere(vec3 pos, vec3 earthCenter, float thickness){
    vec3 posit = (pos - vec3(earthCenter.x, ATMOSPHERE_INNER_RADIUS - EARTH_RADIUS, earthCenter.z))/thickness;

    return posit;
}

// @brief
// Get the point of weather texture to sample

vec2 GetRelativePointToWeatherMap(vec3 positionInDome, vec3 eyePosition, vec3 domeCenter, float radiusDome){
    float r = length(eyePosition - domeCenter);
    float dist = sqrt(radiusDome*radiusDome - r * r);
    vec3 posXmin = eyePosition + vec3(-1.0, 0.0, 0.0) * dist;
    vec3 posZmin = eyePosition + vec3(0.0, 0.0, -1.0) * dist;
    
    posXmin.y = 0.0;
    posZmin.y = 0.0;
    positionInDome.y = 0.0;

    float posx = (positionInDome.x - posXmin.x) / (2.0 * dist); 
    float posz = (positionInDome.z - posZmin.z) / (2.0 * dist);

    //float posx = length(positionInDome - posXmin)/(2 * dist);
    //float posy = length(positionInDome - posZmiz)/(2 * dist);
    vec2 pointSample = vec2(posx, posz);
    return pointSample;
}

/*
 * @brief: HenyeyGreenstein function
 * This function models the atenuation of borders in light scattering
 */
float HenyeyGreenstein(in vec3 inLightVector, in vec3 inViewVector, in float g){
    // Cos of angle
    float cos_angle = dot(normalize(inLightVector), normalize(inViewVector));

    // Define the HenyeyGreenstein function
    return (1.0 - g * g)/(pow(1.0 + g * g - 2.0 * g * cos_angle, 1.50) * 4 * MPI);
}


// ** Define the total energy 
float GetLightEnergy(float density, float probRain, float henyeyGreensteinFactor){
    float beer_laws = 2.0 * exp( -density * probRain);
    float powdered_sugar = 1.0 - exp( -2.0 * density);

    //float henyeyGreensteinFactor = HenyeyGreenstein(lightVector, rayDirection, g);
    float totalEnergy = beer_laws * powdered_sugar * henyeyGreensteinFactor;

    return totalEnergy;
}

vec3 [6] GetNoiseKernel(vec3 lightDirection){
    vec3 maxCompUnitVector;
    if(abs(lightDirection.x) > abs(lightDirection.y) && (abs(lightDirection.x) > abs(lightDirection.z))){
        maxCompUnitVector = vec3(abs(lightDirection.x), 0.0, 0.0);
    }
    else if(abs(lightDirection.y) > abs(lightDirection.x) && (abs(lightDirection.y) > abs(lightDirection.z))){
        maxCompUnitVector = vec3(0.0, abs(lightDirection.y), 0.0);
    }
    else{
        maxCompUnitVector = vec3(0.0, 0.0, abs(lightDirection.z));
    }
    vec3 zcom = cross(lightDirection, maxCompUnitVector);
    vec3 xcom = cross(zcom, lightDirection);
    mat3 sunRotMatrix = mat3(xcom, lightDirection, zcom);
    // Noise light kernell
    vec3 noise_kernel[6] = vec3[6]
    (
        sunRotMatrix * vec3(0.1,    0.25, -0.15),
        sunRotMatrix * vec3(0.2,  0.50,  0.20),
        sunRotMatrix * vec3(-0.2,   0.10, -0.10),
        sunRotMatrix * vec3(-0.05,  0.75,  0.05),
        sunRotMatrix * vec3(-0.1,   1.00,  0.00),
        sunRotMatrix * vec3(0.0,    3.00,  0.00)
    );

    return noise_kernel;
}


// SampleCloudDensity
float SampleCloudDensity(vec3 samplepoint, vec3 weather_data, float relativeHeight, bool ischeap){
    // Add movement to textures
    vec3 wind_direction = vec3(1.0, 0.0, 0.0);
    float cloud_speed = 5.0;

    float cloud_top_offset = 1.0;

    // Skew in wind direction

    samplepoint += relativeHeight * wind_direction * cloud_top_offset * 0.00005;
    samplepoint += (wind_direction + vec3(0.0, 1.0, 0.0))* Time * cloud_speed * 0.001;


    // Init Low Frequency Sampling ...
    vec4 low_frequency_noises = SampleLowFrequencyTexture(samplepoint);

    // Perly Worley noise:
    float low_freq_FBM = low_frequency_noises.y * 0.625 + 
                         low_frequency_noises.z * 0.250 +
                         low_frequency_noises.w * 0.125;

    // Remap with worley noise
    float base_cloud = Remap(low_frequency_noises.x, -(1.0 - low_freq_FBM), 1.0, 0.0, 1.0);

    //Get the gradient density
    float density_height_gradient = GetDensityHeightGradientForPoint(samplepoint, weather_data, relativeHeight);

    // Apply the height function to base cloud
    // Until here the shape of cloud is defined!
    base_cloud = base_cloud * density_height_gradient;

    // Apply the coverage of data
    float cloud_coverage = weather_data.x;
    //base_cloud = clamp(base_cloud, 0.0, 1.0);
    float base_cloud_with_coverage = Remap(base_cloud, cloud_coverage, 1.0, 0.0, 1.0);
    //base_cloud_with_coverage = clamp(base_cloud_with_coverage, 0.0, 1.0);
    // Get more aestheticcal cloud
    base_cloud_with_coverage *= cloud_coverage;

    base_cloud_with_coverage = clamp(base_cloud_with_coverage, 0.0, 1.0);
    
    float final_cloud = base_cloud_with_coverage;

    if(!ischeap && base_cloud_with_coverage > 0.0){
        // sample high frequency texture
        vec3 curl_noise = SampleCurlNoiseTexture(samplepoint.xy);
        samplepoint.xy = samplepoint.xy + curl_noise.xy * (1.0 - relativeHeight);
        vec3 high_frequency_noises = SampleHighFrequencyTexture(samplepoint);

        float high_freq_FBM =     (high_frequency_noises.x * 0.625 )
                                + (high_frequency_noises.y * 0.250 )
                                + (high_frequency_noises.z * 0.125 );

        // TODO: Paper propose other way to compute the height_fraction
        float high_freq_noise_modifier = mix(high_freq_FBM, 1.0 - high_freq_FBM, clamp(relativeHeight * 2.0, 0.0, 1.0));
        
        //high_freq_noise_modifier = clamp(high_freq_noise_modifier, 0.0, 1.0);
        // Erode by remapping:
        final_cloud = Remap(base_cloud_with_coverage, high_freq_noise_modifier * 0.2, 1.0, 0.0, 1.0); 
        final_cloud = clamp(final_cloud, 0.0, 1.0);
    }

    return final_cloud;
    //return clamp(base_cloud, 0.0, 1.0);
    //return base_cloud_with_coverage;
}

/*
 * @brief: SampleCloudAlongCone
 *
 */
float SampleCloudDensityAlongCone(vec3 posInAtm, vec3 rayOrigin, vec3 earthCenter, float stepsize, vec3 noise_kernel[6]){
    float density_along_cone = 0.0;

    float cone_spread_multiplier = stepsize;
    // Lighting ray-march loop
    for(int i = 0; i < 6; i++){
        vec3 lightpos = posInAtm + (cone_spread_multiplier * noise_kernel[i] * float(i));
        vec2 weatherpoint = GetRelativePointToWeatherMap(lightpos, rayOrigin, earthCenter, ATMOSPHERE_OUTER_RADIUS);
        vec3 weather_data = SampleWeatherTexture(weatherpoint);
        
        float relativeHeight = GetRelativeHeightInAtmosphere(lightpos, earthCenter);

        vec3 p = GetPositionInAtmosphere(lightpos, earthCenter, THICK_ATMOSPHERE)*0.5;
        // Sample cloud density:
        density_along_cone += SampleCloudDensity(p, weather_data, relativeHeight, true);
    }

    return density_along_cone;
}


// ** Ray marching algorithm

vec3 RayMarch(vec3 rayOrigin, vec3 startPoint, vec3 endPoint, vec3 rayDirection, vec3 earthCenter, inout float density_inout){
    vec3 colorpixel                 = vec3(0.0);
    float density                   = 0.0;
    float cloud_test                = 0.0;
    int zero_density_sample_count   = 0;
    int sample_cout                 = 128; 
    float thick_                    = length(endPoint - startPoint);
    float stepsize                  = float(thick_/sample_cout);
    float start_                    = length(startPoint - rayOrigin);
    float end_                      = length(endPoint   - rayOrigin);
    vec3 stepSampling               = rayDirection/sample_cout;

    // Henyey Greenstein factor in light Sampling
    // Ray Light direction

    vec3 lightDirection             = normalize(SunLocation - rayOrigin);
    // g: eccentricity 0.2, proposed by paper
    float henyeyGreensteinFactor    = HenyeyGreenstein(lightDirection, rayDirection, 0.2);
    // Get the noise Kernell size 6
    vec3 noise_kernel[6]            = GetNoiseKernel(lightDirection);

    // Start the raymarching loop
    //vec3 samplepoint = GetPositionInAtmosphere(posInAtm, earthCenter, thick_);
    for(float t = start_; t < end_; t += stepsize){
        vec3 posInAtm = rayOrigin + t * rayDirection;
        // Sample the cloud data
        vec3 samplepoint            = GetPositionInAtmosphere(posInAtm, earthCenter, thick_);
        vec2 weatherpoint           = GetRelativePointToWeatherMap(posInAtm, rayOrigin, earthCenter, ATMOSPHERE_OUTER_RADIUS);
        vec3 weather_data           = SampleWeatherTexture(weatherpoint);
        float relativeHeight        = GetRelativeHeightInAtmosphere(posInAtm, earthCenter);
        //float relativeHeight = GetHeightInAtmosphere(posInAtm, earthCenter, startPoint, rayDirection, rayOrigin, thick_);
        // Start with light test 
        if(cloud_test > 0.0){
            //samplepoint = KeepInBox(samplepoint);
            float sampled_density = SampleCloudDensity(samplepoint, weather_data, relativeHeight, false);
            if(sampled_density == 0.0)
                zero_density_sample_count++;
            if(zero_density_sample_count != 6){
                density += sampled_density * 0.05;

                // Start the light sampling
                if(sampled_density != 0.0){
                    // Make the Light Sampling
                    float density_along_light_ray = SampleCloudDensityAlongCone(posInAtm, rayOrigin, earthCenter, stepsize, noise_kernel);
                    float totalEnergy = GetLightEnergy(density_along_light_ray, 0.2, henyeyGreensteinFactor); 
                    float transmitance = 1.0;
                    transmitance = mix(transmitance, totalEnergy, (1.0 - density));

                    colorpixel += vec3(transmitance * 0.5);
                }
                //samplepoint += stepSampling;
                //t += stepsize;
                //posInAtm = rayOrigin +  t * rayDirection;
                //colorpixel += vec3(sampled_density * 0.1);
                //if(sampled_density < 0){
                //    colorpixel = vec3(0.0);
                //    break;
                //}
            }
            else{
                cloud_test = 0.0;
                zero_density_sample_count = 0;
            }
        }
        //Light sampling
        else{
            cloud_test = SampleCloudDensity(samplepoint, weather_data, relativeHeight, true);
            if(cloud_test == 0.0){
                //samplepoint += stepSampling;
                //posInAtm += t * rayDirection;
                //t -= stepsize;
            }       
        }

        if(density >= 1.0){
            density = 1.0;
            break;
        }
        // Define the position in atmosphere 
        //samplepoint = KeepInBox(samplepoint);
    }
    density_inout = density;
    return colorpixel;
}


void main(){
    float aspecRatio    = screenWidth / screenHeight;
    vec3  rayOrigin     = cameraPosition;
    float x             = aspecRatio * (2.0 * gl_FragCoord.x/screenWidth - 1.0);
    float y             = 2.0 * gl_FragCoord.y / screenHeight - 1.0;

    // Get ray Direction

    vec3 rayDirection       = GetRayDirection(cameraFront, cameraRight, cameraUp, x, y);
    
    vec3 earthCenter        = vec3(rayOrigin.x, rayOrigin.y - EARTH_RADIUS, rayOrigin.z);

    vec3 innerIntersection  = GetIntersectionSphereRay(rayOrigin, rayDirection, earthCenter, ATMOSPHERE_INNER_RADIUS);
    vec3 outerIntersection  = GetIntersectionSphereRay(rayOrigin, rayDirection, earthCenter, ATMOSPHERE_OUTER_RADIUS);

    if(dot(vec3(0.0, 1.0, 0.0), rayDirection) < 0.0){
        vec3 colorNearHorizon = vec3(0.54, 0.23, 0.046) * 0.4;
        color =  vec4(colorNearHorizon, 1.0);
        return;
    }

    vec3 skycolor = vec3(0.054687, 0.3, 0.57);
    // ** Perform the color of Sun
    vec3 sunLightDirection = normalize(SunLocation - rayOrigin);
    // ** Dot product to get the intensity of sun at some direction
    float sunFactorEnergy = dot(sunLightDirection, rayDirection);
    sunFactorEnergy = exp(2000.0 * (sunFactorEnergy - 1.0));
    vec3 sunColor = vec3(0.9608, 0.9529, 0.9137);

    skycolor = mix(skycolor, sunColor, sunFactorEnergy);


    float density = 0.0;
    
    vec3 col = RayMarch(rayOrigin, innerIntersection, outerIntersection, rayDirection, earthCenter, density);

    vec3 col_sky_ = mix(skycolor * 1.5, col, density);
    
    //col.x = Remap(skycolor.x, 0.0, col.x, 0.0, 1.0);
    //col.y = Remap(skycolor.y, 0.0, col.y, 0.0, 1.0);
    //col.z = Remap(skycolor.z, 0.0, col.z, 0.0, 1.0);
    //col = skycolor + col;
    //
    //col.x = clamp(col.x, 0.0, 1.0);
    //col.y = clamp(col.y, 0.0, 1.0);
    //col.z = clamp(col.z, 0.0, 1.0);

    
    //vec3 col = SampleCurlNoiseTexture(vec2(x, y)*0.5 + 0.5);
    //vec3 col = vec3(density, density, density);
    //color = vec4(col.xy, 0.0, 1.0);
    color = vec4(col_sky_, 1.0);
}
