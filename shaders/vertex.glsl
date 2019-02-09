#version 330

layout (location = 0) in vec3 pos;

out vec2 screenPos;

uniform float screenWidth;
uniform float screenHeight;
void main(){
    gl_Position = vec4(pos, 1.0);
}