# include <GL/glew.h>
# include <GLFW/glfw3.h>
# include "../inc/Window.hpp"
int main(){
    Window * window = new Window();
    window->Initialize("Window of Clouds");
    glewExperimental = GL_TRUE;
    if(glewInit() != GLEW_OK){
        printf("GLEW initialization fails!!\n");
        delete window;
        return 1;
    }
    glEnable(GL_DEPTH_TEST);
    GLint bufferWidth = window->getBufferWidth();
    GLint bufferHeight = window->getBufferHeight();
    glViewport(0, 0, bufferWidth, bufferHeight);

    while(!window->getShouldClose()){
        glfwPollEvents();
        glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(0);
        window->swapBuffers();
    }
    delete window;
    return 0;
    
}