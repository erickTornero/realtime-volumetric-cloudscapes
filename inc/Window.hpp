#ifndef __WINDOWM__
#define __WINDOWM__

# include <stdio.h>
# include <GL/glew.h>
# include <cstring>
# include <GLFW/glfw3.h>

class Window{
public:
    Window(GLint wWidth = 800, GLint wHeight = 600):width(wWidth), height(wHeight){
        memset(this->keys, false, sizeof(this->keys));
    }
    int Initialize(const char *);

    GLint getBufferWidth(){return this->bufferWidth;}
    GLint getBufferHeight(){return this->bufferHeight;}
    
    bool * getKeys(){return this->keys;}
    GLfloat getXChange();
    GLfloat getYChange();
    bool getShouldClose() {return glfwWindowShouldClose(this->mainWindow);}

    void swapBuffers() {glfwSwapBuffers(this->mainWindow);}

    ~Window();
private:
    GLFWwindow * mainWindow;

    GLint width, height;
    GLint bufferWidth, bufferHeight;

    // ** Handle inputs
    // What key was pressed
    bool keys[1024];

    // ** Handle mouse movement
    GLfloat lastX;
    GLfloat lastY;
    GLfloat xChange;
    GLfloat yChange;

    bool mouseFirstMoved;

    // Callback for handle keys
    void createCallbacks();
    static void handleKeys(GLFWwindow * window, int key, int code, int action, int mode);

    static void handleMouse(GLFWwindow * window, double xPos, double yPos);
};
#endif