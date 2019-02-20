#version 330
// Output color
out vec4 color;
uniform float Time;

// Screen Resolution
uniform float screenWidth;
uniform float screenHeight;
// Camera model
uniform vec3 cameraPosition;
uniform vec3 cameraFront;
uniform vec3 cameraUp;
uniform vec3 cameraRight;
uniform vec3 cameraWorldUp;
uniform mat4 model;

// Parameters
float Epsilon = 0.004;
int MAX_N_STEPS = 200;

float cameraSpeed = 5.0;

float sampleLowFrequency(vec3 point, in vec3 unskewedSamplePoint, in float relativeHeight, in vec3 earthCenter)
{
    //Read in the low-frequency Perlin-Worley noises and Worley noises
    vec4 lowFrequencyNoises = texture(cloudBaseShapeSampler, point);// * 0.8);	// MANIPULATE ME 

    //Build an FBM out of the low-frequency Worley Noises that are used to add detail to the Low-frequency Perlin Worley noise
    float lowFrequencyFBM = (lowFrequencyNoises.g * 0.625) + 
                            (lowFrequencyNoises.b * 0.25)  + 
                            (lowFrequencyNoises.a * 0.125);
    
    // lowFrequencyFBM = clamp(abs(lowFrequencyFBM), 0.0, 1.0);  
    lowFrequencyFBM = clamp(lowFrequencyFBM, 0.0, 1.0);                      

    // Define the base cloud shape by dilating it with the low-frequency FBM
    float baseCloud = remapClamped( lowFrequencyNoises.r, (lowFrequencyFBM - 0.9), 1.0, 0.0, 1.0 );
    
    // TODO: Use weater map for cloud types and blend between them and their densities
    // ------------------ only screws it up; but needed for blending cloud types -----------------------
    // // Get the density-height gradient
    // vec2 weatherSamplePoint = unskewedSamplePoint.xz;// / 50.0;// / 50.0;//50000.0f;
    // vec3 weather_data = texture(weatherMapSampler, weatherSamplePoint).rgb;
    // float cloudType = weather_data.g;
    // float densityHeightGradient = getDensityHeightGradientForPoint(relativeHeight, weather_data.g);

    // // Apply Height function to the base cloud shape
    // baseCloud *= densityHeightGradient * 0.5;// * 0.8;
    // -------------------------------------------------------------------------------------------------

    // Cloud coverage is stored in weather data’s red channel .
    //WHAT EVEN???? --> increasing cloud coverage apparently reduces base_cloud_with_coverage
    float cloud_coverage = 0.6;// weather_data.r;

    // Use remap to apply the cloud coverage attribute.
    float base_cloud_with_coverage = remapClampedBeforeAndAfter ( baseCloud, cloud_coverage, 1.0, 0.0, 1.0);

    // To ensure that the density increases with coverage in an aesthetically pleasing manner
    // Multiply the result by the cloud coverage attribute so that smaller clouds are lighter 
    // and more aesthetically pleasing
    base_cloud_with_coverage *= cloud_coverage;

    return base_cloud_with_coverage;
}


vec3 rayMarch(Ray ray, vec3 earthCenter, in vec3 startPos, in float start_t, in float end_t, in int pixelID, inout float accumDensity)
{
    float _dot = dot(ray.direction, vec3(0.0f, 1.0f, 0.0f));
    float jitterfactor = 1.180f;//4.0f;

    // MANIPULATE ME 
    const float baseDensityFactor = 0.380f;//0.5f; // increase this to get more dense cloud centers    
    const float maxSteps = floor(mix(35.0f, 60.0f, 1.0f - _dot));
	
    const float atmosphereThickness = (end_t - start_t);	
	const float stepSize = (atmosphereThickness / maxSteps);
    float transmittance = 1.0;

    vec3 pos;
    vec3 samplePoint;
    vec3 returnColor = vec3(0.0);
    float baseDensity;

    // Lighting data -------------------------------------------------------------------
    // Henyey-Greenstein
    const vec3 lightDir = normalize(SUN_LOCATION - ray.origin);
    const float cos_angle = dot(normalize(ray.direction), lightDir);
    const float eccentricity = 0.6;
    const float silver_intensity = 0.7;
    const float silver_spread = 0.1;
    const float HG_light = HGModified(cos_angle, eccentricity, silver_intensity, silver_spread);

    // Random unit vectors for your cone sample.
    // These are positioned to be facing the sun 
    // Create random samples within a unit cone facing world up (y direction)
    // Construct TBN matrix towards light
    // Then rotate unit cone towards the sun using TBN

    vec3 maxCompUnitVector;
    if(abs(lightDir[0]) > abs(lightDir[1]) && abs(lightDir[0]) > abs(lightDir[2]))
    {
        maxCompUnitVector = vec3(abs(lightDir[0]), 0.0, 0.0);
    }
    else if(abs(lightDir[1]) > abs(lightDir[0]) && abs(lightDir[1]) > abs(lightDir[2]))
    {
        maxCompUnitVector = vec3(0.0, abs(lightDir[1]), 0.0);
    }
    else
    {
        maxCompUnitVector = vec3(0.0, 0.0, abs(lightDir[2]));
    }
    
    vec3 zComponent = cross(lightDir, maxCompUnitVector);
    vec3 xComponent = cross(zComponent, lightDir);
    mat3 sunRotMatrix = mat3(xComponent, lightDir, zComponent);

    const vec3 noise_kernel[] = 
    {
        sunRotMatrix * vec3(0.1, 0.25, -0.15),
        sunRotMatrix * vec3(0.2, 0.5, 0.2),
        sunRotMatrix * vec3(-0.2, 0.1, -0.1),
        sunRotMatrix * vec3(-0.05, 0.75, 0.05),
        sunRotMatrix * vec3(-0.1, 1.0, 0.0),
        sunRotMatrix * vec3(0.0, 3.0, 0.0),     // One sample should be at distance 3x cone length
    };
    // ---------------------------------------------------------------------------------

	for (float t = start_t; t < end_t; t += stepSize)
	{
		vec3 colorPerSample = vec3(0.0);

        int _index = int(mod((pixelID+int(t)),16.0f));
        vec2 jitterLocation = getJitterOffset(_index, ivec2(75.0f));
		pos = ray.origin + t * (ray.direction + vec3(jitterLocation.x, (jitterLocation.x + jitterLocation.y)*jitterfactor, jitterLocation.y));
        samplePoint = getRelativePositionInAtmosphere(pos, earthCenter);
        samplePoint /= 8.0f; //controls the frequency of how we are sampling the noise texture

		float relativeHeight = getRelativeHeightInAtmosphere(pos, earthCenter, startPos, ray.direction, ray.origin);
        vec3 skewedSamplePoint = skewSamplePointWithWind(samplePoint, relativeHeight);

		baseDensity = sampleLowFrequency(skewedSamplePoint, pos, relativeHeight, earthCenter) * baseDensityFactor; //helps for early termination of rays

		if(baseDensity > 0.0) // Useful to prevent lighting calculations for zero density points
		{
            //Erode Base cloud shape with higher frequency noise (more expensive and so done when we know for sure we are inside the cloud)
            float highFreqDensity = erodeCloudWithHighFrequency(baseDensity*1.4f, ray.direction, skewedSamplePoint, relativeHeight);

            // MANIPULATE ME 
			accumDensity += highFreqDensity * 0.5;

			// Do Lighting calculations with cone sampling
			float densityAlongLight = 0.0;
			int light_samples = 6;

			for(int i = 0; i < light_samples; ++i)
			{
				// Add the current step offset to the sample position
                vec3 lightPos = pos + (stepSize * noise_kernel[i] * float(i));
               	vec3 sampleLightPos =  getRelativePositionInAtmosphere(lightPos, earthCenter);

               	// MANIPULATE ME 
                float currBaseLightDensity = sampleLowFrequency(sampleLightPos, sampleLightPos, relativeHeight, earthCenter);

                if(currBaseLightDensity > 0.0)
                {
                	float currLightDensity = erodeCloudWithHighFrequency(1.5 * currBaseLightDensity, ray.direction, skewedSamplePoint, relativeHeight);
                	densityAlongLight += currLightDensity;
                }
			}

            // ------------------------------------------------------------------------------------------------------------------
            // MANIPULATE ME 
            float brightness = 5.0;
            float totalLightEnergy = GetLightEnergy(relativeHeight, densityAlongLight, baseDensity, HG_light, cos_angle, stepSize, brightness);
            transmittance = mix(transmittance, totalLightEnergy, (1.0 - accumDensity)); 
            colorPerSample = vec3(transmittance);
            
            returnColor += colorPerSample;
		} //end if

		if(accumDensity >= 1.0) 
		{
            accumDensity = 1.0;
            break;
        } //end if

	} //end raymarcher for loop
	return returnColor;
}// end raymarch function

vec3 skewSamplePointWithWind(in vec3 point, inout float height_fraction)
{
    //skew in wind direction
    point += height_fraction * WIND_DIRECTION * CLOUD_TOP_OFFSET * 0.009;
    
    //Animate clouds in wind direction and add a small upward bias to the wind direction
    point += (WIND_DIRECTION + vec3(0.0, 0.1, 0.0)) * CLOUD_SPEED * time.y;
    return point;
}