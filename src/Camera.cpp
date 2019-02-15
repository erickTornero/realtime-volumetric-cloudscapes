# include "../inc/Camera.hpp"

void Camera::update(){
    this->front.x = cos(glm::radians(this->yaw))*cos(glm::radians(this->pitch));
    this->front.y = sin(glm::radians(this->pitch));
    this->front.z = sin(glm::radians(this->yaw))*cos(glm::radians(this->pitch));
    front = glm::normalize(this->front);

    this->right = glm::normalize(glm::cross(this->front, this->worldUp));
    this->up = glm::normalize(glm::cross(this->right, this->front));
}
void Camera::KeyControl(bool * keys, GLfloat dTime){
    GLfloat velocity = movementSpeed*dTime;
    if(keys[GLFW_KEY_W]){
        this->position += this->front*velocity;
    }
    if(keys[GLFW_KEY_S]){
        this->position -= this->front*velocity;
    }
    if(keys[GLFW_KEY_D]){
        this->position += this->right*velocity;
    }
    if(keys[GLFW_KEY_A]){
        this->position -= this->right*velocity;
    }
}

void Camera::MouseControl(GLfloat xChange, GLfloat yChange){
    xChange *= this->turnSpeed;
    yChange *= this->turnSpeed;

    this->yaw += xChange;
    this->pitch += yChange;

    if(this->pitch >= 89.0f)
        this->pitch = 89.0f;
    if(this->pitch < -89.0f)
        this->pitch = -89.0f;
    
    this->update();
}
glm::mat4 Camera::calculateViewMatrix(){
    return glm::lookAt(this->position, this->position + this->front, this->up);
}

void Camera::SetPosition(glm::vec3 newPosition){
    this->position = newPosition;
}
glm::vec3 Camera::getCameraPosition(){
    return this->position;
}
glm::vec3 Camera::getCameraFront(){
    return this->front;
}
glm::vec3 Camera::getCameraUp(){
    return this->up;
}
glm::vec3 Camera::getCameraRight(){
    return this->right;
}
Camera::~Camera(){

}