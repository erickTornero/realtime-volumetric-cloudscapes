# include "../inc/Shader.hpp"

Shader::Shader(){
    this->shaderID = 0;
    this->uniformScreenWidth = 0;
    this->uniformScreenHeight = 0;
    this->uniformCameraPosition = 0;
    //this->uniformModel = 0;
    //this->uniformProjection = 0;
    this->uniformTime = 0;
    //this->uniformMouseXY = 0;
    this->uniformLowFreqTexture = 0;
    this->uniformHighFreqTexture = 0;
    this->uniformWeatherTexture = 0;
    this->uniformGradientStratusTexture = 0;
    this->uniformGradientCumulusTexture = 0;
    this->uniformGradientCumulonimbusTexture = 0;
    this->uniformCurlNoiseTexture = 0;
    this->uniformEarthCenter = 0;
    // Halton vectors
    this->uniformHaltonSeq1 = 0;
    this->uniformHaltonSeq2 = 0;
    this->uniformHaltonSeq3 = 0;
    this->uniformHaltonSeq4 = 0;
}

void Shader::CreateFromString(const char * vertexCode, const char * fragmentCode){
    this->CompileShader(vertexCode, fragmentCode);
}
void Shader::CreateFromFile(const char * fileV, const char * fileF){
    std::string vertexCode = this->ReadFileShader(fileV);
    std::string fragmentCode = this->ReadFileShader(fileF);
    this->CompileShader(vertexCode.c_str(), fragmentCode.c_str());
    printf("\n====================================\n");
    printf("Shader Created Successfully\n");
}

void Shader::CompileShader(const char * vertexCode, const char * fragmentCode){
    this->shaderID = glCreateProgram();
    if(!this->shaderID){
        printf("Error creating shader program!");
        return;
    }
    AddShader(shaderID, vertexCode, GL_VERTEX_SHADER);
    AddShader(shaderID, fragmentCode, GL_FRAGMENT_SHADER);
    GLint result = 0;
    GLchar eLog[1024]={0};

    glLinkProgram(shaderID);
    glGetProgramiv(shaderID, GL_LINK_STATUS, &result);
    if(!result){
        glGetProgramInfoLog(shaderID, sizeof(eLog), NULL, eLog);
        printf("Error Linking program: '%s'\n", eLog);
        return;
    }

    // * Validate program
    glValidateProgram(shaderID);
    glGetProgramiv(shaderID, GL_VALIDATE_STATUS, &result);
    if(!result){
        glGetProgramInfoLog(shaderID, sizeof(eLog), NULL, eLog);
        printf("Error Validating program: '%s'\n", eLog);
        return;
    }
    this->uniformScreenWidth = glGetUniformLocation(shaderID, "screenWidth");
    this->uniformScreenHeight = glGetUniformLocation(shaderID, "screenHeight");
    //uniformXMove = glGetUniformLocation(shader, "xMove");
    //uniformYMove = glGetUniformLocation(shader, "yMove");
    //this->uniformModel = glGetUniformLocation(shaderID, "model");
    //this->uniformProjection = glGetUniformLocation(shaderID, "projection");
    //this->uniformView = glGetUniformLocation(shaderID, "view");
    //this->uniformAmbientColour = glGetUniformLocation(shaderID, "directionalLight.colour");
    //this->uniformAmbientIntensity = glGetUniformLocation(shaderID, "directionalLight.ambientIntensity");
    //this->uniformDirection = glGetUniformLocation(shaderID, "directionalLight.direction");
    //this->uniformDiffuseIntensity = glGetUniformLocation(shaderID, "directionalLight.diffuseIntensity");
    this->uniformCameraPosition = glGetUniformLocation(shaderID, "cameraPosition");
    this->uniformTime = glGetUniformLocation(shaderID, "Time");
    //this->uniformMouseXY = glGetUniformLocation(shaderID, "MouseXY");
    this->uniformCamForward = glGetUniformLocation(shaderID, "cameraFront");
    this->uniformCamUp = glGetUniformLocation(shaderID, "cameraUp");
    this->uniformCamRight = glGetUniformLocation(shaderID, "cameraRight");
    this->uniformLowFreqTexture = glGetUniformLocation(shaderID, "lowFrequencyTexture");
    this->uniformHighFreqTexture = glGetUniformLocation(shaderID, "highFrequencyTexture");
    this->uniformWeatherTexture = glGetUniformLocation(shaderID, "WeatherTexture");
    //this->uniformGradientStratusTexture = glGetUniformLocation(shaderID, "GradientStratusTexture");
    //this->uniformGradientCumulusTexture = glGetUniformLocation(shaderID, "GradientCumulusTexture");
    //this->uniformGradientCumulonimbusTexture = glGetUniformLocation(shaderID, "GradientCumulonimbusTexture");
    this->uniformCurlNoiseTexture = glGetUniformLocation(shaderID, "CurlNoiseTexture");
    this->uniformEarthCenter = glGetUniformLocation(shaderID, "EarthCenter");

    this->uniformHaltonSeq1 = glGetUniformLocation(shaderID, "HaltonSequence1");
    this->uniformHaltonSeq2 = glGetUniformLocation(shaderID, "HaltonSequence2");
    this->uniformHaltonSeq3 = glGetUniformLocation(shaderID, "HaltonSequence3");
    this->uniformHaltonSeq4 = glGetUniformLocation(shaderID, "HaltonSequence4");
}

void Shader::AddShader(GLuint theProgram, const char * shaderCode, GLenum shaderType){
    GLuint theShader = glCreateShader(shaderType);

    const GLchar * theCode[1];
    theCode[0] = shaderCode;
    
    GLint codeLength[1];
    codeLength[0] = strlen(shaderCode);

    glShaderSource(theShader, 1, theCode, codeLength);
    glCompileShader(theShader);

    GLint result = 0;
    GLchar eLog[1024]={0};

    glGetShaderiv(theShader, GL_COMPILE_STATUS, &result);
    if(!result){
        glGetShaderInfoLog(theShader, sizeof(eLog), NULL, eLog);
        printf("Error Compiling the %d shader: '%s'\n", shaderType, eLog);
        return;
    }
    glAttachShader(theProgram, theShader);
}

//GLint Shader::GetModelLocation(){
//    return this->uniformModel;
//}
//GLint Shader::GetProjectionLocation(){
//    return this->uniformProjection;
//}
//GLint Shader::GetViewLocation(){
//    return this->uniformView;
//}/*
//GLint Shader::GetAmbientColourLocation(){
//    return this->uniformAmbientColour;
//}
//GLint Shader::GetIntensityLocation(){
//    return this->uniformAmbientIntensity;
//}
//GLint Shader::GetDiffuseIntensityLocation(){
//    return this->uniformDiffuseIntensity;
//}
//GLint Shader::GetDirectionLocation(){
//    return this->uniformDirection;
//}*/
GLint Shader::GetScreenWidthLocation(){
    return this->uniformScreenWidth;
}
GLint Shader::GetScreenHeightLocation(){
    return this->uniformScreenHeight;
}

GLint Shader::GetCameraPositionLocation(){
    return this->uniformCameraPosition;
}
GLint Shader::GetTimeLocation(){
    return this->uniformTime;
}
//GLint Shader::GetMouseXYLocation(){
//    return this->uniformMouseXY;
//}
GLint Shader::GetCamForwardLocation(){
    return this->uniformCamForward;
}
GLint Shader::GetCamUpLocation(){
    return this->uniformCamUp;
}
GLint Shader::GetCamRightLocation(){
    return this->uniformCamRight;
}
GLint Shader::GetLowFreqTextureLocation(){
    return this->uniformLowFreqTexture;
}
GLint Shader::GetWeatherTextureLocation(){
    return this->uniformWeatherTexture;
}
GLint Shader::GetGradientStratusTextureLocation(){
    return this->uniformGradientStratusTexture;
}
GLint Shader::GetGradientCumulusTextureLocation(){
    return this->uniformGradientCumulusTexture;
}
GLint Shader::GetGradientCumulonimbusTextureLocation(){
    return this->uniformGradientCumulonimbusTexture;
}
GLint Shader::GetHighFreqTextureLocation(){
    return this->uniformHighFreqTexture;
}
GLint Shader::GetCurlNoiseTextureLocation(){
    return this->uniformCurlNoiseTexture;
}
GLint Shader::GetEarthCenterLocation(){
    return this->uniformEarthCenter;
}
// Halton vectors:
GLint Shader::GetHaltonSeq1Location(){
    return this->uniformHaltonSeq1;
}
GLint Shader::GetHaltonSeq2Location(){
    return this->uniformHaltonSeq2;
}
GLint Shader::GetHaltonSeq3Location(){
    return this->uniformHaltonSeq3;
}
GLint Shader::GetHaltonSeq4Location(){
    return this->uniformHaltonSeq4;
}
void Shader::UseShader(){
    glUseProgram(this->shaderID);
}
void Shader::ClearShader(){
    if(this->shaderID != 0){
        glDeleteProgram(this->shaderID);
        this->shaderID = 0;
    }
    this->uniformScreenWidth = 0;
    this->uniformScreenHeight = 0;
    //this->uniformModel = 0;
    //this->uniformProjection = 0; 
    this->uniformCameraPosition = 0;
    this->uniformTime = 0;
    //this->uniformMouseXY = 0;

    this->uniformCamForward = 0;
    this->uniformCamUp = 0;
    this->uniformCamRight = 0;
    this->uniformLowFreqTexture = 0;
    this->uniformWeatherTexture = 0;
    this->uniformGradientStratusTexture = 0;
    this->uniformGradientCumulusTexture = 0;
    this->uniformGradientCumulonimbusTexture = 0;
    this->uniformHighFreqTexture = 0;
    this->uniformCurlNoiseTexture = 0;
    this->uniformEarthCenter = 0;
    
    //Haton
    this->uniformHaltonSeq1 = 0;
    this->uniformHaltonSeq2 = 0;
    this->uniformHaltonSeq3 = 0;
    this->uniformHaltonSeq4 = 0;
}
std::string Shader::ReadFileShader(const char * fileLoc){
    std::string content;
    std::ifstream file(fileLoc, std::ios::in);
    if(!file.is_open()){
        printf("Fail Opening File ...\n");
        return "";
    }
    std::string line = "";
    while(!file.eof()){
        std::getline(file, line);
        content.append(line + "\n");
    }
    file.close();
    return content;
}
Shader::~Shader(){
    this->ClearShader();
}