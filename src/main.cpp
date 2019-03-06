
# include <GL/glew.h>
# include <GLFW/glfw3.h>
# include <glm/glm.hpp>
# include <glm/gtc/matrix_transform.hpp>
# include <glm/gtc/type_ptr.hpp>
# include <iostream>
# include "../inc/Window.hpp"
# include "../inc/Shader.hpp"
# include "../inc/Mesh.hpp"
# include "../inc/Camera.hpp"
# include "../inc/Texture.hpp"

GLfloat deltaTime = 0.0f;
GLfloat lastTime  = 0.0f;
void calcAverageNormals(unsigned int * indexes, unsigned int indiceCount, GLfloat * vertices, unsigned int verticeCount, unsigned int vLength, unsigned int normalOffset){
    for(size_t i = 0; i < indiceCount; i+=3){
        // ** Getting the indexes of Vertexes of a triangle
        unsigned int in0 = indexes[i] * vLength;
        unsigned int in1 = indexes[i + 1] * vLength;
        unsigned int in2 = indexes[i + 2] * vLength;
        // ** Getting the normal vector of that surface (triangle)
        // ** First Get 2 vectors that define tha plane of the triangle
        glm::vec3 v1(vertices[in1] - vertices[in0], vertices[in1 + 1] - vertices[in0 + 1], vertices[in1 + 2] - vertices[in0 + 2]);
        glm::vec3 v2(vertices[in2] - vertices[in0], vertices[in2 + 1] - vertices[in0 + 1], vertices[in2 + 2] - vertices[in0 + 2]);
        // ** Then the cross product return the normal vector to the surface.
        glm::vec3 normal = glm::cross(v1, v2);
        // ** Normalize this vector as norm equals to 1
        normal = glm::normalize(normal);

        // ** Go to the position of normal components by take in account the normall Offset (in this case 5)
        in0 += normalOffset; in1 += normalOffset; in2 += normalOffset;
        // ** Then Update the normal components of each Surface.
        vertices[in0] += normal.x; vertices[in0 + 1] += normal.y; vertices[in0 + 2] += normal.z;
        vertices[in1] += normal.x; vertices[in1 + 1] += normal.y; vertices[in1 + 2] += normal.z;
        vertices[in2] += normal.x; vertices[in2 + 1] += normal.y; vertices[in2 + 2] += normal.z;
    }

    for(size_t i = 0; i < verticeCount/vLength; i++){
        unsigned int nOffset = i * vLength + normalOffset;
        glm::vec3 vec(vertices[nOffset], vertices[nOffset + 1], vertices[nOffset + 2]);
        vec = glm::normalize(vec);
        vertices[nOffset] = vec.x;
        vertices[nOffset + 1] = vec.y;
        vertices[nOffset + 2] = vec.z;
    }
}

Mesh * CreateTriangle(){
    // Create a set of indexes, that specifies how to take each of the four triangles
    unsigned int indexes[] = {
        0, 3, 1,    //Define the first triangle
        1, 3, 2,
        0, 2, 3,
        0, 1, 2
    };
    // Four points can define a piramid in the 3D Space
    // So we define 4 vertex.
    // ** u, v define the coordinates of the texture
    // ** from (0,0) to (1,1)
    GLfloat vertices[]={
    //**  x     y       z       u   v       nx      ny      nz
        -1.0f, -1.0f, 0.0f,    0.0f, 0.0f,  0.0f,   0.0f,   0.0f,
        0.0f, -1.0f, 1.0f,     0.5f, 0.0f,  0.0f,   0.0f,   0.0f,
        1.0f, -1.0f, 0.0f,     1.0f, 0.0f,  0.0f,   0.0f,   0.0f,
        0.0f, 1.0f, 0.0f,       0.5f, 1.0f, 0.0f,   0.0f,   0.0f
    };
    calcAverageNormals(indexes, 12, vertices, 32, 8, 5);
    Mesh *m = new Mesh();
    m->CreateMesh(vertices, indexes, 32, 12);
    return m;
    //objects.push_back(m);
}
Mesh * CreateQuad(){
    unsigned int indexes[] = {
        0, 1, 2,
        0, 1, 3
    };
    GLfloat vertices[] = {
        // x    y   z
        -1.0, -1.0, -1.0,
         1.0,  1.0, -1.0,
        -1.0,  1.0, -1.0,
         1.0, -1.0, -1.0
    };
    Mesh * mesh = new Mesh();
    mesh->CreateMesh(vertices, indexes, 12, 6);
    return mesh;
}

int main(){
    Window * window = new Window(640, 480);
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
    std::cout<<"W> "<<bufferWidth<<std::endl;
    std::cout<<"H> "<<bufferHeight<<std::endl;
    glViewport(0, 0, bufferWidth, bufferHeight);
    Mesh * quad = CreateQuad();
    //Mesh * quadstatic = CreateQuad();
    Shader * shader = new Shader();
    Texture * lowfreqTexture = new Texture("textures/LowFrequency3DTexture.tga");
    Texture * highFreqTexture = new Texture("textures/HighFrequency3DTexture.tga");
    Texture * weatherTexture = new Texture("textures/weather_t1.tga");
    //Texture * gradientStratus = new Texture("textures/gradient_stratus2d.tga");
    //Texture * gradientCumulus = new Texture("textures/gradient_cumulus2d.tga");
    //Texture * gradientCumulonimbus = new Texture("textures/gradient_cumulonimbus2d.tga");
    Texture * curlNoiseTexture = new Texture("textures/curlNoise.png");
    // Load 3d Texture in RGBA format
    lowfreqTexture->LoadTexture3D();
    // Load 3d texture
    highFreqTexture->LoadTexture3D();
    // Load 2D texture in RGB format
    weatherTexture->LoadTexture();
    // Load 1D Textures in Grayscale for Height gradient functions
    //gradientStratus->LoadTexture2DGray();
    //gradientCumulus->LoadTexture2DGray();
    //gradientCumulonimbus->LoadTexture2DGray();
    curlNoiseTexture->LoadTextureA();

    //std::cout<<"width> "<<lowfreqTexture->
    shader->CreateFromFile("shaders/vertex.glsl", "shaders/RayMarching2.glsl");
    Camera * camera = new Camera(glm::vec3(0.0, 0.0, -2.0), glm::vec3(0.0, 1.0, 0.0), -90.0, 0.0, 500.0, 0.03);
    // Seth the initial position of earth
    glm::vec3 earthCenter = camera->getCameraPosition() - glm::vec3(0.0, 6378000.0, 0.0);
    //glm::mat4 projection = glm::perspective(45.0f, (GLfloat)bufferWidth/(GLfloat)bufferHeight, 0.1f, 100.0f);
    GLint indexTexture = 0;
    while(!window->getShouldClose()){
        GLfloat now = glfwGetTime();
        deltaTime = now - lastTime;
        lastTime = now;
        glfwPollEvents();
        camera->KeyControl(window->getKeys(), deltaTime);
        camera->MouseControl(window->getXChange(), window->getYChange());
        glClearColor(0.74, 0.84, 0.02, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        shader->UseShader();


        glUniform1f(shader->GetScreenWidthLocation(), GLfloat(bufferWidth));
        glUniform1f(shader->GetScreenHeightLocation(), GLfloat(bufferHeight));
        /*
         *   Camera plane model
        */
        //glm::mat4 model(1.0f);
        //model = glm::translate(model, glm::vec3(-1.2f, -2.0, -2.5f));
        //glm::vec4 postemp = model*glm::vec4(0.0,0.0,0.0, 1.0);
        //std::cout<<"("<<camera->getCameraPosition().x<<", "<<camera->getCameraPosition().y<<", "<<camera->getCameraPosition().z<<")"<<std::endl;
        //camera->SetPosition(glm::vec3(postemp.x, postemp.y, postemp.z) - glm::vec3(0.0, 0.0, -2.0));
        //glUniformMatrix4fv(shader->GetModelLocation(), 1, GL_FALSE, glm::value_ptr(model));
        //glUniformMatrix4fv(shader->GetProjectionLocation(), 1, GL_FALSE, glm::value_ptr(projection));
        //glUniformMatrix4fv(shader->GetViewLocation(), 1, GL_FALSE, glm::value_ptr(camera->calculateViewMatrix()));
        glUniform3fv(shader->GetCameraPositionLocation(), 1, glm::value_ptr(camera->getCameraPosition()));
        glUniform3fv(shader->GetCamForwardLocation(),1, glm::value_ptr(camera->getCameraFront()));
        glUniform3fv(shader->GetCamUpLocation(),1, glm::value_ptr(camera->getCameraUp()));
        glUniform3fv(shader->GetCamRightLocation(),1, glm::value_ptr(camera->getCameraRight()));
        //glUniform1f()
        glUniform1f(shader->GetTimeLocation(), now);
        // Pass the earth Center as variable:
        glUniform3fv(shader->GetEarthCenterLocation(), 1, glm::value_ptr(earthCenter));
        //glUniform2fv(shader->GetMouseXYLocation(), 1, glm::value_ptr(glm::vec2(window->getXChange(), window->getYChange())));
        quad->RenderMesh();
        lowfreqTexture->UseTexture3D(shader->GetLowFreqTextureLocation(), indexTexture++);
        highFreqTexture->UseTexture3D(shader->GetHighFreqTextureLocation(), indexTexture++);
        //glUniform1i(shader->GetLowFreqTextureLocation(), 0);
        //glActiveTexture(GL_TEXTURE0);
        //glBindTexture(GL_TEXTURE_3D, lowfreqTexture->GetID());
        
        weatherTexture->UseTexture(shader->GetWeatherTextureLocation(), indexTexture++);
        //glUniform1i(shader->GetWeatherTextureLocation(), 1);
        //glActiveTexture(GL_TEXTURE1);
        //glBindTexture(GL_TEXTURE_2D, weatherTexture->GetID());
        
        //gradientStratus->UseTexture(shader->GetGradientStratusTextureLocation(), indexTexture++);
        //gradientCumulus->UseTexture(shader->GetGradientCumulusTextureLocation(), indexTexture++);
        //gradientCumulonimbus->UseTexture(shader->GetGradientCumulonimbusTextureLocation(), indexTexture++);
        curlNoiseTexture->UseTexture(shader->GetCurlNoiseTextureLocation(), indexTexture++);
        std::cout<<camera->getCameraPosition().x<<", "<<camera->getCameraPosition().y<<", "<<camera->getCameraPosition().z<<std::endl;
        
        
        //glUniform1i(shader->GetGradientCumulonimbusTextureLocation(), 2);
        //glActiveTexture(GL_TEXTURE2);
        //glBindTexture(GL_TEXTURE_2D, gradientCumulonimbus->GetID());
        
        //glUniform1i(shader->GetGradientStratusTextureLocation(), 3);
        //glActiveTexture(GL_TEXTURE3);
        //glBindTexture(GL_TEXTURE_1D, gradientStratus->GetID());
        /*
         *  Static Quad model
         */ 
        /*model = glm::mat4(1.0f);
        model = glm::translate(model, glm::vec3(0.0, -5.0f, -4.0f));
        model = glm::rotate(model, (glm::mediump_float)90*0.0175f, glm::vec3(0.8f, 0.0f, 0.0f));
        glUniformMatrix4fv(shader->GetModelLocation(), 1, GL_FALSE, glm::value_ptr(model));
        quadstatic->RenderMesh();
        */
        //glUniform()
        indexTexture = 0;
        
        glUseProgram(0);
        window->swapBuffers();
        
        
    }
    
    printf("Releasing Resources\n");
    delete window;
    delete shader;
    delete camera;
    delete quad;
    delete lowfreqTexture;
    delete weatherTexture;
    delete highFreqTexture;
    //delete gradientStratus;
    //delete gradientCumulus;
    //delete gradientCumulonimbus;

    return 0;
    
}
