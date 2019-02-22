#ifndef __SHADERM__
#define __SHADERM__

# include <stdio.h>
# include <string.h>
# include <iostream>
# include <fstream>

#include <GL/glew.h>

class Shader{
public:
    Shader();

    void CreateFromString(const char * vertexCode, const char * fragmentCode);
    void CreateFromFile(const char * fileV, const char * fileF);
    GLint GetProjectionLocation();
    GLint GetModelLocation();
    GLint GetViewLocation();
    //GLint GetIntensityLocation();
    //GLint GetAmbientColourLocation();
    //GLint GetDiffuseIntensityLocation();
    //GLint GetDirectionLocation();
    GLint GetScreenWidthLocation();
    GLint GetScreenHeightLocation();

    GLint GetCameraPositionLocation();
    GLint GetTimeLocation();
    GLint GetMouseXYLocation();
    GLint GetCamForwardLocation();
    GLint GetCamUpLocation();
    GLint GetCamRightLocation();

    // 3D texture sampler
    GLint GetLowFreqTextureLocation();
    // Weather Texture 2D sampler
    GLint GetWeatherTextureLocation();
    // Density height function 1D texture 1D;
    GLint GetGradientStratusTextureLocation();
    GLint GetGradientCumulusTextureLocation();
    GLint GetGradientCumulonimbusTextureLocation();

    void UseShader();
    
    ~Shader();
private:
    GLuint shaderID;
    GLint uniformProjection, uniformModel, uniformView;
    GLint uniformCameraPosition;
    GLint uniformTime;
    GLint uniformMouseXY;
    GLint uniformCamForward;
    GLint uniformCamUp;
    GLint uniformCamRight;
    // 3D Texture 128x128x128 in RGBA format  
    GLint uniformLowFreqTexture;
    // 2D Texture 512x512 in RGB format
    GLint uniformWeatherTexture;
    // 1D Texture 300 x 1 Grayscale format for height density functions
    GLint uniformGradientStratusTexture;
    GLint uniformGradientCumulusTexture;
    GLint uniformGradientCumulonimbusTexture;
    //GLint uniformAmbientIntensity, uniformAmbientColour, uniformDiffuseIntensity, uniformDirection;
    GLint uniformScreenWidth, uniformScreenHeight;
    void CompileShader(const char * vertexCode, const char * fragmentCode);
    void AddShader(GLuint theProgram, const char * shaderCode, GLenum shaderType);
    void ClearShader();
    std::string ReadFileShader(const char * fileLoc);
};
#endif