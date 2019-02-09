#version 330

layout (location = 0) in vec3 pos;

out vec2 screenPos;

uniform float screenWidth;
uniform float screenHeight;
void main(){
    screenPos.x = (pos.x + 1)*(screenWidth  * 2.);
    screenPos.y = (pos.y + 1)*(screenHeight * 2.);
    gl_Position = vec4(pos, 1.0);
}