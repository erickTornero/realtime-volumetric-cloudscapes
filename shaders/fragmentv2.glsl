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