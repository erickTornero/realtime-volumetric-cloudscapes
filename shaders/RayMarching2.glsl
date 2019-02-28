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

// Gradient heigh density

float GetGradientHeightFactor(float height, int cloudtype){
    float timewidthup, starttimeup, starttimedown;
    if(cloudtype == 0){
        timewidthup = 0.08;
        starttimeup = 0.08;
        starttimedown = 0.2;
    }
    else if(cloudtype == 1){
        timewidthup = 0.14;
        starttimeup = 0.1;
        starttimedown = 0.5;
    }
    else if(cloudtype == 2){
        timewidthup = 0.2;
        starttimeup = 0.10;
        starttimedown = 0.7;
    }

    float factor = 2 * MPI/(2 * timewidthup);

    float density = 0.0;
    // Gradient functions dependent on the cloudtype
    if(height < starttimeup)
        density = 0.0;
    else if(height < starttimeup + timewidthup)
        density = 0.5 * sin(factor * height - MPI/2.0 - factor * starttimeup) + 0.5;
    else if(height < starttimedown)
        density = 1.0;
    else if(height < starttimedown + timewidthup)
        density = 0.5 * sin(factor * height - MPI/2.0 - factor * (starttimedown + timewidthup)) + 0.5;
    else
        density = 0.0;

    return density;
}

float GetDensityHeightGradientForPoint(vec3 point, vec3 weather_data){
    float cloudt = weather_data.z;
    // Cloud type: {0: Stratus, 1: Cumulus, 2: Cumulonimbus}
    int cloudtype = 1;
    if(cloudt < 0.1)
        cloudtype = 0;
    else if(cloudt > 0.9)
        cloudtype = 2;
    
    // Relative Height from [0 - 1]
    float relativeHeight = point.y;
    //float relativeHeight = GetRelativeHeightInAtmosphere(point, earthCenter);

    // Get gradient function defined for three type of clouds
    return GetGradientHeightFactorighFactor(relativeHeight, cloudtype);
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

// SampleCloudDensity
float SampleCloudDensity(vec3 samplepoint, vec3 weather_data, bool ischeap){
    vec4 low_frequency_noises = SampleLowFrequencyTexture(samplepoint);

    // Perly Worley noise:
    float low_freq_FBM = low_frequency_noises.y * 0.625 + 
                         low_frequency_noises.z * 0.250 +
                         low_frequency_noises.w * 0.125;

    // Remap with worley noise
    float base_cloud = Remap(low_frequency_noises.x, -(1.0 - low_freq_FBM), 1.0, 0.0, 1.0);

    //Get the gradient density
    float density_height_gradient = GetDensityHeightGradientForPoint(samplepoint, weather_data);

    // Apply the height function to base cloud
    // Until here the shape of cloud is defined!
    base_cloud *= density_height_gradient;

    // Apply the coverage of data
    //float cloud_coverage = weather_data.x;
    //float base_cloud_with_coverage = Remap(base_cloud, cloud_coverage, 1.0, 0.0, 1.0);

    // Get more aestheticcal cloud
    //base_cloud_with_coverage *= cloud_coverage;

    //float final_cloud = base_cloud_with_coverage;

    //if(!ischeap){
        // sample high frequency texture

    //}

    return base_cloud;

}
// ** Ray marching algorithm

float RayMarch(vec3 startPoint, vec3 rayDirection, vec3 rayDirection){
    float density = 0.0;
    float cloud_test = 0.0;
    int zero_density_sample_count = 0;

    vec3 stepSampling = rayDirection/sample_cout;
    // Start the raymarching loop
    for(int i = 0; i < sample_cout; i++){
        // Start with light test 
        if(cloud_test > 0.0){
            float sampled_density = SampleCloudDensity(samplepoint, weather_data, false);
            if(sampled_density = 0.0)
                zero_density_sample_count++;
            if(zero_density_sample_count != 6){
                density += sampled_density;
                samplepoint += stepSampling;
            }
            else{
                cloud_test = 0.0;
                zero_density_sample_count = 0;
            }
        }
        // Light sampling
        else{
            cloud_test = SampleCloudDensity(samplepoint, weather_data, true);
            if(cloud_test == 0.0){
                samplepoint += stepSampling;
            }       
        }

        samplepoint = KeepInBox(samplepoint);
    }
    return density;
}


void main(){
    
}