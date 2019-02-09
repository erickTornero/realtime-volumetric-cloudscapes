# include <GL/glew.h>
# include <GLFW/glfw3.h>
# include <glm/glm.hpp>
# include "../inc/Window.hpp"
# include "../inc/Shader.hpp"
# include "../inc/Mesh.hpp"

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
    Mesh * quad = CreateQuad();
    Shader * shader = new Shader();
    shader->CreateFromFile("shaders/vertex.glsl", "shaders/fragment.glsl");
    while(!window->getShouldClose()){
        glfwPollEvents();
        glClearColor(0.95, 0.95, 0.95, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        shader->UseShader();

        quad->RenderMesh();
        glUseProgram(0);
        window->swapBuffers();
    }
    delete window;
    return 0;
    
}