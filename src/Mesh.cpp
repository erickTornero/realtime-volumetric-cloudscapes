#include "../inc/Mesh.hpp"

Mesh::Mesh(){
    this->VAO = 0;
    this->VBO = 0;
    this->IBO = 0;
    this->indexCount = 0;
}
void Mesh::CreateMesh(GLfloat * vertices, unsigned int * indexes, unsigned int numVertex, unsigned int numIndexes){
    this->indexCount = numIndexes;

    glGenVertexArrays(1, &this->VAO);
    glBindVertexArray(this->VAO);

    glGenBuffers(1, &this->IBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->IBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indexes[0])*numIndexes, indexes, GL_STATIC_DRAW);
     
    glGenBuffers(1, &this->VBO);
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices[0])*numVertex, vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(0);
    //glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(vertices[0])*8, (void*)(sizeof(vertices[0])*3));
    //glEnableVertexAttribArray(1);
    //glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(vertices[0])*8, (void*)(sizeof(vertices[0])*5));
    //glEnableVertexAttribArray(2);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}
void Mesh::RenderMesh(){
    glBindVertexArray(this->VAO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->IBO);
    glDrawElements(GL_TRIANGLES, this->indexCount, GL_UNSIGNED_INT, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

}

void Mesh::ClearMesh(){
    if(this->IBO != 0){
        glDeleteBuffers(1, &this->IBO);
        this->IBO = 0;
    }
    if(this->VBO != 0){
        glDeleteBuffers(1, &this->VBO);
        this->VBO = 0;
    }
    if(this->VAO != 0){
        glDeleteVertexArrays(1, &this->VAO);
        this->VAO = 0;
    }
    this->indexCount = 0;
}
Mesh::~Mesh(){
    this->ClearMesh();
}
