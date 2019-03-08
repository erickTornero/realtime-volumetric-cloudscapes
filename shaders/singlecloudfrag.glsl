#version 330

out vec4 color;

uniform vec3 cameraPosition;
uniform vec3 cameraFront;
uniform vec3 cameraUp;
uniform vec3 cameraRight;
// Earth center position
uniform vec3 EarthCenter;

uniform float screenWidth;
uniform float screenHeight;
uniform float Time;

uniform sampler3D lowFrequencyTexture;
uniform sampler3D highFrequencyTexture;
uniform sampler2D WeatherTexture;

uniform sampler2D CurlNoiseTexture;


const float MPI = 3.14159265;

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
    vec3 ray = 2 * front +  right * x + up * y;
    return normalize(ray);
}
// ** Remap function defined by the author of paper
float Remap(float original_value, float original_min, float original_max, float new_min, float new_max){
    return new_min + (((original_value - original_min)/(original_max - original_min)) * (new_max - new_min));
}

bool ComputeInterSection(vec3 center, float radius, vec3 rayOrigin, vec3 rayDirection, inout vec3 stp, inout vec3 ep){
    vec3 distV = rayOrigin - center;
    float dist = length(distV);

    float a = length(rayDirection);
    float b = 2 * dot(distV, rayDirection);
    float c = dist * dist - radius * radius;

    float discriminat = b * b - 4 * a * c;
    float t1, t2;
    if(discriminat < 0) {
        return false;
    }
    else{
        t1 = (-b + sqrt(discriminat))/(2.0*a);
        t2 = (-b - sqrt(discriminat))/(2.0*a);
        if(t1 < 0.0)
            t1 = 0.0;
        if(t2 < 0.0)
            t2 = 0.0;

        if(t2 < t1){
            float tmp = t1;
            t1 = t2;
            t2 = tmp;
        }
        stp = rayOrigin + t1 * rayDirection;
        ep  = rayOrigin + t2 * rayDirection;
        return true;
    }
}

// SampleCloudDensity
float SampleCloudDensity(vec3 samplepoint, bool ischeap){
    // Add movement to textures
    //vec3 wind_direction = vec3(1.0, 0.0, 0.0);
    //float cloud_speed = 10.0;
//
    //float cloud_top_offset = 5.0;
//
    //// Skew in wind direction
//
    //samplepoint += relativeHeight * wind_direction * cloud_top_offset * 0.005;
    //samplepoint += (wind_direction + vec3(0.0, 1.0, 0.0))* Time * cloud_speed * 0.001;


    // Init Low Frequency Sampling ...
    vec4 low_frequency_noises = SampleLowFrequencyTexture(samplepoint);

    // Perly Worley noise:
    float low_freq_FBM = low_frequency_noises.g * 0.625 + 
                         low_frequency_noises.b * 0.250 +
                         low_frequency_noises.a * 0.125;


    //low_freq_FBM = clamp(low_freq_FBM, 0.0, 1.0);


    // Remap with worley noise
    float base_cloud = Remap(low_frequency_noises.r, -(1.0 - low_freq_FBM), 1.0, 0.0, 1.0);
    //base_cloud = clamp(base_cloud, 0.0, 1.0);
    //Get the gradient density
    //float density_height_gradient = GetDensityHeightGradientForPoint(samplepoint, weather_data, relativeHeight);

    // Apply the height function to base cloud
    // Until here the shape of cloud is defined!
    //base_cloud = base_cloud * density_height_gradient;

    // Apply the coverage of data
    //float cloud_coverage = weather_data.r;
    //cloud_coverage = clamp(cloud_coverage, 0.0, base_cloud)
    //base_cloud = clamp(base_cloud, cloud_coverage, 1.0);
    //float base_cloud_with_coverage = Remap(base_cloud, cloud_coverage, 1.0, 0.0, 1.0);
    //base_cloud_with_coverage = clamp(base_cloud_with_coverage, 0.0, 1.0);
    //// Get more aestheticcal cloud
    //base_cloud_with_coverage *= cloud_coverage;

    //base_cloud_with_coverage = clamp(base_cloud_with_coverage, 0.0, 1.0);
    
    //float final_cloud = base_cloud_with_coverage;

    float final_cloud = base_cloud;
    if(!ischeap ){
        // sample high frequency texture
        vec3 curl_noise = SampleCurlNoiseTexture(samplepoint.xy);
        samplepoint.xy = samplepoint.xy + curl_noise.xy * (1.0 - 0.5);
        vec3 high_frequency_noises = SampleHighFrequencyTexture(samplepoint * 0.1);

        float high_freq_FBM =     (high_frequency_noises.r * 0.625 )
                                + (high_frequency_noises.g * 0.250 )
                                + (high_frequency_noises.b * 0.125 );

        high_freq_FBM = clamp(high_freq_FBM, 0.0, 1.0);
        // TODO: Paper propose other way to compute the height_fraction
        //float high_freq_noise_modifier = mix(high_freq_FBM, 1.0 - high_freq_FBM, clamp(relativeHeight * 10.0, 0.0, 1.0));
        //high_freq_noise_modifier *= 0.35 * exp(-cloud_coverage * 0.75);
        //base_cloud_with_coverage = clamp(base_cloud_with_coverage, 0.2*high_freq_noise_modifier, 1.0);
        //high_freq_noise_modifier = clamp(high_freq_noise_modifier, 0.0, 1.0);
        //base_cloud_with_coverage = clamp(base_cloud_with_coverage, 0.2 * high_freq_noise_modifier, 1.0);
        // Erode by remapping:
        final_cloud = Remap(final_cloud, 0.2, 1.0, 0.0, 1.0);
        //final_cloud = clamp(final_cloud, 0.0, 1.0);
    }

    return final_cloud;
    //return clamp(base_cloud, 0.0, 1.0);
    //return base_cloud_with_coverage;
}


vec3 RayMarch(vec3 rayOrigin, vec3 startPoint, vec3 endPoint, vec3 rayDirection, inout float density_inout){
    vec3 colorpixel                 = vec3(0.0);
    float density                   = 0.0;
    float cloud_test                = 0.0;
    int zero_density_sample_count   = 0;
    //int sample_cout                 = 128;
    int sample_cout                 = 128;// - 64 * int(dot(rayDirection, vec3(0.0, 1.0, 0.0))); 
    float thick_                    = length(endPoint - startPoint);
    float stepsize                  = float(thick_/sample_cout);
    float start_                    = length(startPoint - rayOrigin);
    float end_                      = length(endPoint   - rayOrigin);
    if(end_ < start_){
        float tmp = start_;
        start_ = end_;
        end_ = tmp;
    }
    vec3 stepSampling               = rayDirection/float(sample_cout);

    // Henyey Greenstein factor in light Sampling
    // Ray Light direction
    vec3 testsize                   = rayDirection/float(sample_cout);

    //vec3 lightDirection             = normalize(SunLocation - rayOrigin);
    // g: eccentricity 0.2, proposed by paper
    //float henyeyGreensteinFactor    = HenyeyGreenstein(lightDirection, rayDirection, 0.6);
    // Get the noise Kernell size 6
    //vec3 noise_kernel[6]            = GetNoiseKernel(lightDirection);

    vec3 samplepoint = vec3(0.5, 0.5, 0.0);
    // Start the raymarching loop
    //vec3 samplepoint = GetPositionInAtmosphere(posInAtm, earthCenter, thick_);
    for(float t = start_; t < end_; t += stepsize){
        //vec2 jitterLoc = getJitterOffset(0, vec2(80.0));
        vec3 posInAtm = rayOrigin + t * (rayDirection);// + vec3(jitterLoc.x, (jitterLoc.x + jitterLoc.y)*4.0, jitterLoc.y));
        // Sample the cloud data
        //vec3 samplepoint            = GetPositionInAtmosphere(posInAtm, earthCenter, thick_);
        //vec2 weatherpoint           = GetRelativePointToWeatherMap(posInAtm, rayOrigin, earthCenter, ATMOSPHERE_OUTER_RADIUS);
        //vec3 weather_data           = SampleWeatherTexture(weatherpoint);
        //float relativeHeight        = GetRelativeHeightInAtmosphere(posInAtm, earthCenter);
        //float relativeHeight = GetHeightInAtmosphere(posInAtm, earthCenter, startPoint, rayDirection, rayOrigin, thick_);
        // Start with light test 
        if(cloud_test > 0.0){
            samplepoint += testsize;
            samplepoint /=1.0;
            //samplepoint = KeepInBox(samplepoint);
            float sampled_density = SampleCloudDensity(samplepoint * 200.0, false);
            if(sampled_density == 0.0)
                zero_density_sample_count++;
            if(zero_density_sample_count != 6){
                density += sampled_density ;
                colorpixel += vec3(sampled_density * 0.1);
                // Start the light sampling
                /*
                if(sampled_density != 0.0){
                    // Make the Light Sampling
                    float precipitation;
                    if(weather_data.z > 0.9) precipitation = 5.0;
                    else if(weather_data.z < 0.1) precipitation = 0.05;
                    else precipitation = 2.0;

                    float density_along_light_ray = SampleCloudDensityAlongCone(posInAtm, rayOrigin, earthCenter, stepsize, noise_kernel);
                    float totalEnergy = GetLightEnergy(density_along_light_ray, precipitation, henyeyGreensteinFactor); 
                    float transmitance = 1.0;
                    transmitance = mix(density, totalEnergy, (1.0 - density_along_light_ray));
                    //transmitance = clamp(transmitance, 0.0, 1.0);
                    colorpixel += vec3(transmitance * 0.4);
                    //colorpixel += vec3(sampled_density);
                    //colorpixel += vec3(totalEnergy * 1.0);
                }*/
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
            cloud_test = SampleCloudDensity(samplepoint, true);
            if(cloud_test == 0.0){
                //samplepoint += stepSampling;
                //posInAtm += t * rayDirection;
                t -= stepsize;
            }       
        }

        if(density >= 1.0){
            density = 1.0;
            break;
        }
        // Define the position in atmosphere 
        //samplepoint = KeepInBox(samplepoint);
    }
    //colorpixel = vec3(density);
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
    
    //vec3 earthCenter        = vec3(rayOrigin.x, rayOrigin.y - EARTH_RADIUS, rayOrigin.z);
    //vec3 earthCenter        = EarthCenter;
    //vec3 innerIntersection  = GetIntersectionSphereRay(rayOrigin, rayDirection, earthCenter, ATMOSPHERE_INNER_RADIUS);
    //vec3 outerIntersection  = GetIntersectionSphereRay(rayOrigin, rayDirection, earthCenter, ATMOSPHERE_OUTER_RADIUS);

    /*
    if(dot(vec3(0.0, 1.0, 0.0), rayDirection) < 0.0){
        vec3 colorNearHorizon = vec3(0.54, 0.23, 0.046) * 0.4;
        color =  vec4(colorNearHorizon, 1.0);
        return;
    }
    else if(dot(vec3(0.0, 1.0, 0.0), rayDirection) < 0.06){
        vec3 skycolor = getAtmosphereColorPhysical(rayDirection, SunLocation - rayOrigin, SunIntensity);
        skycolor *= 0.620;
        //vec3 colorNearHorizon = vec3(0.1, 0.1, 0.1);
        //color = vec4(colorNearHorizon, 1.0);
        color = vec4(skycolor, 1.0);
        return;
    }
    vec3 skycolor = getAtmosphereColorPhysical(rayDirection, SunLocation - rayOrigin, SunIntensity);
    skycolor *= 0.62;
    //vec3 skycolor = vec3(0.054687, 0.3, 0.57);
    // ** Perform the color of Sun
    vec3 sunLightDirection = normalize(SunLocation - rayOrigin);
    // ** Dot product to get the intensity of sun at some direction
    float sunFactorEnergy = dot(sunLightDirection, rayDirection);
    sunFactorEnergy = exp(2000.0 * (sunFactorEnergy - 1.0));
    vec3 sunColor = vec3(0.9608, 0.9529, 0.9137);
    */

    //skycolor = mix(skycolor, 2.0 * sunColor, sunFactorEnergy);

    
    float density = 0.0;
    vec3 spherecenter = vec3(0.0, 0.0, -40.0);
    float radius = 10.0;
    vec3 startPoint, endPoint;

    bool thereis = ComputeInterSection(spherecenter, radius, rayOrigin, rayDirection, startPoint, endPoint);
    
    vec3 col;
    if(thereis){
        col = RayMarch(rayOrigin, startPoint, endPoint, rayDirection, density);
        //density *= smoothstep(0.0, 1.0, min(1.0, Remap(rayDirection.y, 0.06, 0.4, 0.0, 1.0)));
        //col = vec3(density);
    }
    else
        col = vec3(0.0);
    
    //vec3 sk_c = mix(col,  10.0 * sunColor, sunFactorEnergy);
    //vec3 col_sky_ = mix(skycolor, col, density);
    
    //vec3 col_sky_ = mix(skycolor * 1.5, col, density);
    
    //col.x = Remap(skycolor.x, 0.0, col.x, 0.0, 1.0);
    //col.y = Remap(skycolor.y, 0.0, col.y, 0.0, 1.0);
    //col.z = Remap(skycolor.z, 0.0, col.z, 0.0, 1.0);
    //col = skycolor + col;
    //
    //col.x = clamp(col.x, 0.0, 1.0);
    //col.y = clamp(col.y, 0.0, 1.0);
    //col.z = clamp(col.z, 0.0, 1.0);

    //vec4 cc = SampleLowFrequencyTexture(vec3(vec2(x,y), 0.0));
    //color = vec4(cc.w, cc.w, cc.w, 1.0);
    //vec3 col = SampleCurlNoiseTexture(vec2(x, y)*0.5 + 0.5);
    //vec3 col = vec3(density, density, density);
    //color = vec4(col.xy, 0.0, 1.0);

    // Debug WeatherMap
    //vec2 wp = GetRelativePointToWeatherMap(outerIntersection, rayOrigin, earthCenter, ATMOSPHERE_OUTER_RADIUS);


    //vec3 c = SampleWeatherTexture(wp);
    //vec3 c = SampleWeatherTexture(vec2(x, y) * 0.5 + 0.5);
    //color = vec4(c, 1.0);
    //color = vec4(col_sky_, 1.0);

    color = vec4(col, 1.0);
}
