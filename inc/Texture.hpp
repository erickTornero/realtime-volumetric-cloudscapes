# ifndef __TEXTUREM__
# define __TEXTUREM__

# include <GL/glew.h>

#include "stb_image.h"

class Texture{
public:
    Texture(const char * fileLoc = ""):textureID(0),width(0), height(0), bitDepht(0), fileLocation(fileLoc){};

    bool LoadTexture();
    // ** Work with alpha channel
    bool LoadTextureA();
    // ** 3D Texture
    bool LoadTexture3D();
    // ** 1D Texture Load in Grayscale
    bool LoadTexture1D();

    void UseTexture(GLint, GLint);
    void UseTexture1D(GLint, GLint);
    void UseTexture3D(GLint, GLint);

    void ClearTexture();
    ~Texture();

private:
    GLuint textureID;
    int width, height, bitDepht;

    const char * fileLocation;
};

# endif