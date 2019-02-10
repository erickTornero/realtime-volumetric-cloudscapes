#version 330

layout (location = 0) in vec3 pos;

out vec2 screenPos;
//out vec3 cameraPosition;
uniform float screenWidth;
uniform float screenHeight;
uniform vec3 cameraPosition;

uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
void main(){
    // screenPos.x = (pos.x + 1)*(screenWidth  * 2.);
    // screenPos.y = (pos.y + 1)*(screenHeight * 2.);
    //cameraPosition = (model * vec4(pos, 1.0)).xyz;
    gl_Position = vec4(pos, 1.0);
    //gl_Position = vec4(cameraPosition + vec3(0, 0, -10.0), 1.0);
    //gl_Position = projection * view * model * vec4(pos, 1.0);
}