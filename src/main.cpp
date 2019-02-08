# include <GL/glew.h>
# include <GLFW/glfw3.h>
# include "../inc/Window.hpp"
# include "../inc/Shader.hpp"
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

    Shader * shader = new Shader();
    shader->CreateFromFile("shaders/vertex.glsl", "shaders/fragment.glsl");
    while(!window->getShouldClose()){
        glfwPollEvents();
        glClearColor(0.95, 0.95, 0.95, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        shader->UseShader();
        glUseProgram(0);
        window->swapBuffers();
    }
    delete window;
    return 0;
    
}