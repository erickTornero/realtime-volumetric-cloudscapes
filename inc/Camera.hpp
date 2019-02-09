#ifndef __CAMERAM__
#define __CAMERAM__

# include <GL/glew.h>

# include <glm/glm.hpp>
# include <glm/gtc/matrix_transform.hpp>
# include <GLFW/glfw3.h>

class Camera{
public:
    Camera(glm::vec3 position, glm::vec3 up, GLfloat yaw, GLfloat pitch, GLfloat msp, GLfloat tsp):position(position), worldUp(up), yaw(yaw), pitch(pitch), front(glm::vec3(0.0f, 0.0f, -1.0f)), movementSpeed(msp), turnSpeed(tsp){
        this->update();
    }

    glm::mat4 calculateViewMatrix();

    void KeyControl(bool * keys, GLfloat dTime);
    void MouseControl(GLfloat xChange, GLfloat yChange);
    ~Camera();

private:
    glm::vec3 position;
    glm::vec3 front;
    glm::vec3 up;
    glm::vec3 right;
    glm::vec3 worldUp;

    GLfloat yaw;
    GLfloat pitch;

    GLfloat movementSpeed;
    GLfloat turnSpeed;

    void update();
};


#endif