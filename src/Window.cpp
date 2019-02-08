# include "../inc/Window.hpp"


int Window::Initialize(const char * nameW = "My Window"){
    //Initialize GLFW
    if(!glfwInit()){
        printf("GLFW Initializing fails\n");
        glfwTerminate();
        return 1;
    }
    // Setup GLFW window properties
    // OpenGL Version
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    // Gives an error if lower version of OpenGL is used
    // Core profile = No backward compatibility
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Allow forward Compability
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    
    //Create a windows
    this->mainWindow = glfwCreateWindow(this->width, this->height, nameW, NULL, NULL);
    // Get Buffer size information
    glfwGetFramebufferSize(this->mainWindow, &this->bufferWidth, &this->bufferHeight);

    // Set the context for GLEW to use
    glfwMakeContextCurrent(this->mainWindow);

    //handle key & mouse
    this->createCallbacks();
    glfwSetInputMode(this->mainWindow, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    //test
    glfwSetWindowUserPointer(this->mainWindow, this);
}

GLfloat Window::getXChange(){
    GLfloat theChange = this->xChange;
    this->xChange = 0.0f;
    return theChange;
}

GLfloat Window::getYChange(){
    GLfloat theChange = this->yChange;
    this->yChange = 0.0f;
    return theChange;
}

void Window::handleKeys(GLFWwindow * window, int key, int code, int action, int mode){
    Window * theWindow = static_cast<Window*>(glfwGetWindowUserPointer(window));
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
    if(key >= 0 && key < 1024){
        if(action == GLFW_PRESS){
            theWindow->keys[key] = true;
        }
        else if(action == GLFW_RELEASE){
            theWindow->keys[key] = false;
        }
    }
}
void Window::handleMouse(GLFWwindow * window, double xPos, double yPos){
    Window * theWindow = static_cast<Window*>(glfwGetWindowUserPointer(window));

    if(theWindow->mouseFirstMoved){
        theWindow->lastX = xPos;
        theWindow->lastY = yPos;
        theWindow->mouseFirstMoved = false;
    }

    theWindow->xChange = xPos - theWindow->lastX;
    theWindow->yChange = theWindow->lastY - yPos;

    theWindow->lastX = xPos;
    theWindow->lastY = yPos;

    // sprintf("x:%.6f, y:%.6f\n", theWindow->xChange, theWindow->yChange);
}
void Window::createCallbacks(){
    glfwSetKeyCallback(this->mainWindow, handleKeys);
    glfwSetCursorPosCallback(this->mainWindow, handleMouse);
}
Window::~Window(){
    if(this->mainWindow != nullptr){
        glfwDestroyWindow(mainWindow);
        glfwTerminate();
    }
}
