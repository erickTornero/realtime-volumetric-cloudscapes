# define STB_IMAGE_IMPLEMENTATION
# include "../inc/Texture.hpp"


bool Texture::LoadTexture(){
    unsigned char * texData = stbi_load(this->fileLocation, &this->width, &this->height, &this->bitDepht, 0);
    if(!texData){
        printf("Failed Loading the file texture\n");
        return false;
    }
    glGenTextures(1, &this->textureID);
    glBindTexture(GL_TEXTURE_2D ,this->textureID);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, this->width, this->height, 0, GL_RGB, GL_UNSIGNED_BYTE, texData);
    glGenerateMipmap(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, 0);

    stbi_image_free(texData);
    printf("\n==================================\n");
    printf("Texture 2D in RGB format Loaded!\n");
    printf("width>%i\n", this->width);
    printf("height>%i\n", this->height);
    printf("depth>%i\n", this->bitDepht);
    
    return true;
}

bool Texture::LoadTextureA(){
    unsigned char * texData = stbi_load(this->fileLocation, &this->width, &this->height, &this->bitDepht, 0);
    if(!texData){
        printf("Failed Loading the file texture\n");
        return false;
    }
    glGenTextures(1, &this->textureID);
    glBindTexture(GL_TEXTURE_2D ,this->textureID);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, this->width, this->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
    glGenerateMipmap(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, 0);

    stbi_image_free(texData);

    printf("\n==================================\n");
    printf("Texture 2D in RGBA format Loaded!\n");
    printf("width>%i\n", this->width);
    printf("height>%i\n", this->height);
    printf("depth>%i\n", this->bitDepht);

    return true;
}

bool Texture::LoadTexture3D(){
    unsigned char * texData3D = stbi_load(this->fileLocation, &this->width, &this->height, &this->bitDepht, 0);
    if(!texData3D){
        printf("Failed Loading the file texture\n");
        return false;
    }
    glGenTextures(1, &this->textureID);
    glBindTexture(GL_TEXTURE_3D, this->textureID);
    //
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, 128, 128, 128, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData3D);
    stbi_image_free(texData3D);

    printf("\n==================================\n");
    printf("Texture 3D in RGBA format Loaded!\n");
    printf("width>%i\n", this->width);
    printf("height>%i\n", this->height);
    printf("depth>%i\n", this->bitDepht);
    return true;
}
void Texture::UseTexture3D(GLint textureLocation, GLint indexTexture){

    glActiveTexture(GL_TEXTURE0 + indexTexture);
    glUniform1i(textureLocation, indexTexture);
    glBindTexture(GL_TEXTURE_3D, this->textureID);
}


void Texture::UseTexture(GLint textureLocation, GLint indexTexture){
    glActiveTexture(GL_TEXTURE0 + indexTexture);
    glUniform1i(textureLocation, indexTexture);
    glBindTexture(GL_TEXTURE_2D, this->textureID);
}
void Texture::ClearTexture(){
    glDeleteTextures(1, &textureID);
    this->textureID = 0;
    this->width = 0;
    this->height = 0;
    this->bitDepht = 0;
    this->fileLocation = "";
}

Texture::~Texture(){
    this->ClearTexture();
}