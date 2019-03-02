#version 330

out vec4 color;

uniform vec3 cameraPosition;
uniform vec3 cameraFront;
uniform vec3 cameraUp;
uniform vec3 cameraRight;

uniform float screenWidth;
uniform float screenHeight;
uniform float Time;
// Variable for sample textures
uniform sampler3D lowFrequencyTexture;
uniform sampler3D highFrequencyTexture;
uniform sampler2D WeatherTexture;

uniform sampler2D GradientCumulusTexture;
uniform sampler2D GradientCumulonimbusTexture;
uniform sampler2D GradientStratusTexture;
uniform sampler2D CurlNoiseTexture;

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

const float MPI = 3.14159265;

// ****************** Light Variables ***************

const vec3 sun_location = vec3(0.0, ATMOSPHERE_OUTER_RADIUS*0.9, - ATMOSPHERE_OUTER_RADIUS *0.9);

// Remap from book

float Remap(float original_value, float original_min, float original_max, float new_min, float new_max){
    return new_min + (((original_value - original_min)/(original_max - original_min)) * (new_max - new_min));
}

float densityHeightFactor(float height, int typecloud){
    float pi = 3.14159265;
    float timewidthup, starttimeup, starttimedown;
    if(typecloud == 0){
        timewidthup = 0.08;
        starttimeup = 0.08;
        starttimedown = 0.2;
    }
    else if(typecloud == 1){
        timewidthup = 0.14;
        starttimeup = 0.1;
        starttimedown = 0.5;
    }
    else if(typecloud == 2){
        timewidthup = 0.2;
        starttimeup = 0.10;
        starttimedown = 0.7;
    }

    float factor = 2*pi/(2*timewidthup);
    
    float density = 0.0;
    if(height < starttimeup)
        density = 0.0;
    else if(height < starttimeup + timewidthup)
        density = 0.5*sin(factor * height - pi/2.0 - factor * starttimeup) + 0.5;
    else if(height < starttimedown)
        density = 1.0;
    else if(height < starttimedown + timewidthup)
        density =  0.5 * sin(factor * height - pi/2.0 - factor * (starttimedown + timewidthup)) + 0.5;
    else
        density = 0.0;

    return density;
    //for i in range(100):
	//step = i/100.0
	//if step < starttimeup:
	//	y.append(0.0)
	//elif step < starttimeup + timewidthup:
	//	y.append(0.5 * math.sin(factor * step - math.pi/2 - 39.27*starttimeup) + 0.5)
	//elif step < starttimedown:
	//	y.append(1.0)
	//elif step < starttimedown + timewidthup:
	//	y.append(0.5 * math.sin(factor * step - math.pi/2 - factor * (starttimedown + timewidthup)) + 0.5)
	//else:
	//	y.append(0.0)
}
/* @brief
 * Density height function based on wheather map
 */
vec3 SampleWeatherTexture(vec2 point){
    vec4 weatherdata = texture(WeatherTexture, point);
    return weatherdata.xyz;
}
vec3 SampleCurlNoiseTexture(vec2 point){
    vec4 curlnoisedata = texture(CurlNoiseTexture, point);
    return curlnoisedata.xyz;
}

vec3 SampleHighFrequencyTexture(vec3 point){
    vec4 highFreqTexture = texture(highFrequencyTexture, point);
    return highFreqTexture.xyz;
}
float GetDensityHeightGradientForPoint(vec3 weatherdata, float height){
    // Sample the texture pn point
    //vec4 weatherdata = texture(WeatherTexture, point);
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

    float heightFactor = densityHeightFactor(height, cloudtype);
    //vec4 heightFactor = vec4(1.0);
    //if(cloudtype == 0)
    //    heightFactor = texture(GradientStratusTexture, 1.0 - vec2(height, height));
    //else if(cloudtype == 1)
	//    heightFactor = texture(GradientCumulusTexture, 1.0 - vec2(height, height));
    //else
	//   heightFactor = texture(GradientCumulonimbusTexture, 1.0 - vec2(height, height));
    
    // TODO: for moment just cloud coverage is been using.
    return heightFactor;

}

/*
 * @brief
 * Get the texture weather map intesection
 * by mapping the texture over all the dome
 */
vec2 GetPointToSampleInDome(vec3 positionInDome, vec3 eyePosition, vec3 domeCenter, float radiusDome){
    float r = length(eyePosition - domeCenter);
    float dist = sqrt(radiusDome * radiusDome - r * r);
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
/*
float remapClamped(float value, float old_min, float old_max, float new_min, float new_max){
    float v = new_min +((value - old_min)/ (old_max - old_min))*(new_max - new_min);
    return clamp(v, new_min, new_max);
}
float remapClampPrevPost(float value, float old_min, float old_max, float new_min, float new_max){
    value = clamp(value, old_min, old_max);
    return remapClamped(value, old_min, old_max, new_min, new_max);
}*/

vec3 GetRayDirection(vec3 front, vec3 right, vec3 up, float x, float y){
    vec3 ray = 5 * front + right * x + up * y;
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

    return clamp(posInAtm/atmThick, 0.0, 1.0);
}

// Simple way to get relative height
// Get Relative cloud in thick atmosphere
float GetRelativeHeightInAtmosphere(vec3 position, vec3 earthCenter){
    float distCenter2Point = length(position - earthCenter);
    float relativeDistance = distCenter2Point - ATMOSPHERE_INNER_RADIUS;

    return relativeDistance/ATMOSPHERE_THICKNESS;
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

float sampleCloudDensity(vec3 p, vec2 weather_point, float relativeHeight, bool ischeap){
    // Simulate movement of clouds
    vec3 wind_direction = vec3(1.0, 0.0, 0.0);
    float cloud_speed = 10.0;
    float cloud_top_offset = 1.0;
    p += relativeHeight * wind_direction * cloud_top_offset * 0.0001;
    p += (wind_direction + vec3(0.0, 1.0, 0.0)) * Time * cloud_speed * 0.01;
    //weather_point += wind_direction.xy*0.0003 * Time;
    // Read the low frequency Perlin-Worley & Worley noise
    vec4 low_frequency_noises = texture(lowFrequencyTexture, p);
    float low_freq_FBM = low_frequency_noises.y * 0.625 + 
                        low_frequency_noises.z * 0.250 +
                        low_frequency_noises.w * 0.125;
    float base_cloud = Remap(low_frequency_noises.x, -(1.0 - low_freq_FBM), 1.0, 0.0, 1.0);

    // TODO: Base on whether texture
    vec3 weatherdata = SampleWeatherTexture(weather_point);
    float density_height_gradient = GetDensityHeightGradientForPoint(weatherdata, relativeHeight);

    base_cloud = base_cloud * density_height_gradient;
    float cloud_coverage = weatherdata.x;
    // Experimental clamping before & after remaping
    base_cloud = clamp(base_cloud, 0.0, 1.0);
    float base_cloud_coverage = Remap(base_cloud, cloud_coverage, 1.0, 0.0, 1.0);
    base_cloud_coverage = clamp(base_cloud_coverage, 0.0, 1.0);
    //base_cloud_coverage = clamp(base_cloud, 0.0, 1.0);
    base_cloud_coverage = base_cloud_coverage * cloud_coverage;

    float final_cloud = base_cloud_coverage;
    if(!ischeap){
        //if(base_cloud_coverage > 0.0){
        // high frequency working
        vec3 curl_noise = SampleCurlNoiseTexture(p.xy);
        p.xy += curl_noise.xy * (1.0 - relativeHeight);
        vec3 high_frequency_noises = SampleHighFrequencyTexture(p * 0.1);
        float high_freq_FBM =(high_frequency_noises.x * 0.625 + high_frequency_noises.y * 0.25 + high_frequency_noises.z * 0.125, 0.0, 1.0);
        float high_freq_noise_modifier = mix(high_freq_FBM, 1.0 - high_freq_FBM, clamp(10.0 * relativeHeight, 0.0, 1.0));

        final_cloud = Remap(base_cloud_coverage, 0.02 * high_freq_noise_modifier, 1.0, 0.0, 1.0);
        final_cloud = clamp(final_cloud, 0.0, 1.0);
        //}
    }
    
    return final_cloud;
    //return clamp(base_cloud_coverage, 0.0, 1.0);
}
float sampleCloudDensityAlongCone(vec3 position, float stepsize, vec3 noise_kernel[6], vec3 earthCenter, float atmosphereThickness){
    float cone_spread_multiplier = stepsize;
    float density_along_cone = 0.0;
    for(int i = 0; i< 6; i++){
        // Att the current step offset to the sample position
        vec3 lightPos = position + (cone_spread_multiplier * noise_kernel[i] * float(i));
        float rHeight = GetRelativeHeightInAtmosphere(position, earthCenter);
        vec2 weatherpoint = GetPointToSampleInDome(lightPos, earthCenter + vec3(0.0, EARTH_RADIUS, 0.0), earthCenter, ATMOSPHERE_OUTER_RADIUS);
        //vec3 weatherdata = SampleWeatherTexture(weatherpoint);
        vec3 p = GetPositionInAtmosphere(position, earthCenter, atmosphereThickness);
        //if(density_along_view_cone < 0.3){
            density_along_cone += sampleCloudDensity(p, weatherpoint, rHeight, false);
        //}
        //else{
         //   density_along_cone += sampleCloudDensity(p, weatherdata, relativeHeight, true);
        //}
    }

    return density_along_cone;
}
/*
 * @brief: Return the energy transmited along the ray:
 * d: density along the ray
 * p: absortion multiplier for rain
 * g eccentricity in light direction
 * theta: angle between view & light rays
 */
float GetLightEnergy(float d, float p, float g, float theta){
    float lightEnergy = 1.0;
    float beerlaw = exp(-d * p);
    float powder_sugar_effect = 1.0 - exp(-2.0 * d);

    lightEnergy = beerlaw * powder_sugar_effect;
    // Henyey greensteing phase function
    float greenstein = 2.0 * (1 - g * g)/(4 * MPI * (1 + g * g - 2 * g * sqrt(cos(theta)*cos(theta)*cos(theta))));

    lightEnergy *= greenstein;
    return lightEnergy;
}

/*float sampleLowFrequencyTexture(vec3 pointSample){
    // 3D texture sampled 4 channels RGBA
    vec4 sampledTexture = texture(lowFrequencyTexture, pointSample);
    float valueFBM = sampledTexture.y * 0.625 + sampledTexture.z * 0.25 + sampledTexture.w * 0.125;

    valueFBM = clamp(valueFBM, 0.0, 1.0);
    float coverage = 0.75;
    float baseCloud = remapClamped(sampledTexture.x, (valueFBM - 0.9), 1.0, 0.0, 1.0);
    float baseCloudCoverage = remapClampPrevPost(baseCloud, coverage, 1.0, 0.0, 1.0);

    baseCloudCoverage *= coverage;

    return baseCloudCoverage;
}*/

vec3 RayMarching(vec3 rayOrigin, vec3 rayDirection, vec3 innerIntersection, vec3 earthCenter, float start_atm, float end_atm){
    float atmosphereThickness = end_atm - start_atm;
    // max N Steps in accordance to the author of paper
    float maxNSteps = 128.0;
    float incrementAtm = atmosphereThickness/maxNSteps;
    
    // ** Lighting variables
    vec3 lightDirection = normalize(rayOrigin - sun_location);
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
    
    float cloud_test = 0.0;
    int zero_density_sample_count = 0;
    float acummDensity = 0.0;
    vec3 colorPixel = vec3(0.0);
    // Start Ray marching
    for(float t = start_atm; t < end_atm; t += incrementAtm){
        vec3 colorSample = vec3(0.0);

        // Compute position to sample in atmosphere
        vec3 pos = rayOrigin + t * rayDirection;
        //positionInDome, vec3 eyePosition, vec3 domeCenter, float radiusDome
        vec2 weatherdata = GetPointToSampleInDome(pos, rayOrigin, earthCenter, ATMOSPHERE_OUTER_RADIUS);
        float relativeHeight = GetHeightInAtmosphere(pos, earthCenter, innerIntersection, rayDirection, rayOrigin, atmosphereThickness);

        // Warn: Hardcoded variable
        vec3 pointSampled = GetPositionInAtmosphere(pos, earthCenter, atmosphereThickness);//vec3(0.6, 0.0, 0.6);
        pointSampled /=4.0;
        //vec2 weatherPoint = GetPointToSampleInDome(pos, rayOrigin, earthCenter, ATMOSPHERE_OUTER_RADIUS);
        // no cheap sampling
        if(cloud_test > 0.0){
            float sampled_density = sampleCloudDensity(pointSampled, weatherdata, relativeHeight, false);
            if(sampled_density == 0.0)
                zero_density_sample_count++;
            if(zero_density_sample_count != 6){
                acummDensity += sampled_density;
                // Start with Light sampling
                
                float densityAlongLight = 0.0;
                densityAlongLight = sampleCloudDensityAlongCone(pos, incrementAtm, noise_kernel, earthCenter, atmosphereThickness);
                int light_samples = 6;
                //for(int i = 0; i < light_samples; i++){
                //    // Add the current step offset
                //    vec3 lightpos = posInAtm + (incrementAtm * noise_kerel[i] * float(i))
                //    vec3 sampleLightPos = GetPositionInAtmosphere(lightpos, earthCenter, atmosphereThickness);
//
                //    float currBaseLightDensity = sampleCloudDensity(sampleLightPos, sampleLightPos.xy, relativeHeight, earthCenter, true);
                //    if(currBaseLightDensity > 0.0){
                //        float curLightDensity = sampleCloudDensity(currBaseLightDensity, weatherdata, relativeHeight, earthCenter, false);
                //        densityAlongLight += curLightDensity;
                //    }
                //}
                float brightness = 5.0;
                float theta = dot(normalize(lightDirection), rayDirection);
                float totalLightEnergy = GetLightEnergy(densityAlongLight, 0.8, 0.2, theta);
                float transmitance = mix(1.0, totalLightEnergy, (1.0 - acummDensity));
                /// End ligth treatment
                //colorPixel += vec3(sampled_density * 0.1);
                colorPixel += vec3(transmitance*0.03);
            }
            else{
                cloud_test = 0.0;
                zero_density_sample_count = 0;
            }
        }// Cheap sampling
        else{
            cloud_test = sampleCloudDensity(pointSampled, weatherdata, relativeHeight, true);
            if(cloud_test > 0.0){
                t -= incrementAtm;
            }
        }
        //float baseCloud = sampleCloudDensity(pointSampled, weatherPoint, relativeHeight, true);
        
        //float baseCloud = sampleLowFrequencyTexture(pointSampled);
        //float baseCloud = 0.5;
        
        //if(baseCloud > 0.0){
        //    acummDensity += baseCloud *0.006;
        //    //acummDensity += baseCloud * 0.1;
        //    colorPixel += vec3(baseCloud * 0.05);
        //    //break;
        //}
        if(acummDensity >= 1.0){
            acummDensity = 1.0;
            break;
        }
        //if(baseCloud < 0){
        //    colorPixel = vec3(0.0);
        //    break;
        //}

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
    float finalLength   = length(outerIntersection - rayOrigin);

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
    //vec4 weatherdata = texture(WeatherTexture, vec2(x,y)/2.0 + 0.5);
    //vec4 gradient = texture(GradientCumulonimbusTexture, -vec2(y,y)/2.0 + 0.5);
    //vec4 gradient = texture(GradientStratusTexture, y/2.0 + 0.5);
    
    //vec3 col = vec3(0, 0, weatherdata.x);
    //vec3 col = SampleHighFrequencyTexture(vec3(x, y, 0.1)).xyz;
    
    //vec3 col = SampleCurlNoiseTexture(vec2(x, y)/2.0 + 0.5);

    //float typecloud = weatherdata.z;
    //float red = 0.0;
    //float green = 0.0;
    //float blue = 0.0;
    //if(weatherdata.z < 0.55 && weatherdata.z > 0.45)
    //    green = 1.0;
    //else if(weatherdata.z > 0.9)
    //    blue = 1.0;
    //else if(weatherdata.z < 0.1)
    //   red = 1.0;
    //
    //vec3 col = vec3(red, green, blue);
    
    
    vec3 col = RayMarching(rayOrigin, rayDirection, innerIntersection, earthCenter, initialLength, finalLength);
    //col = vec3(1.0 - col.x, 1.0 - col.y, 1.0 - col.z);
    //float modcolor = length(col);
    //if(modcolor < 0.4)
    //    col = vec3(0.5273, 0.808, 0.922);
    color = vec4(col, 1.0);
}
