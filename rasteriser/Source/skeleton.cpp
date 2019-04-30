#include <iostream>
#include <time.h>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include "glm/ext.hpp"
#include <omp.h>
#include <map>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::vec2;
using glm::ivec2;

#define SCREEN_WIDTH 512
#define SCREEN_HEIGHT 512
#define FULLSCREEN_MODE false

vec4 cameraPos(0, 0, -2.001, 1); //removing the 0.001 will cause a crash to occour
vec4 cameraRot(0, 0, 0, 1);
vec4 cameraDir(0, 0, 1, 0);

vec4 lightPos(0,-0.7,0,1);
float theta = 0.0f;
//vec4 lightPos(-0.8, -0.5, -0.8, 1.0);
vec3 lightPower = 7.0f*vec3(1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1);
vec3 globalReflectance(1.5);
vec4 currentNormal;

vec3 band1(1.5);
vec3 band2(3);

struct Pixel {
    int x;
    int y;
    float z;
    vec4 pos3d;
};

struct Vertex {
    vec4 position;
};

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
float focalLength = SCREEN_WIDTH/2;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

/*
//For debugging
void DrawVertecies(screen* screen, vector<Vertex> vertices){
    vec3 color(0.75f, 0.15f, 0.15f);

    int vs = vertices.size();
    vector<Pixel> vertexPixels(vs);
    for(int i = 0; i < vs; i++) {
        VertexShader(vertices[i], vertexPixels[i]);
        //cout << "x: " << vertexPixels[i].x << "\n" << "y: " << vertexPixels[i].y << "\n";
    }

    for (uint32_t i = 0; i < vertices.size(); i++){
      SafePutPixelSDL(screen, vertexPixels[i].x, vertexPixels[i].y-1, color);
      SafePutPixelSDL(screen, vertexPixels[i].x, vertexPixels[i].y+1, color);
      SafePutPixelSDL(screen, vertexPixels[i].x, vertexPixels[i].y, color);
      SafePutPixelSDL(screen, vertexPixels[i].x+1, vertexPixels[i].y-1, color);
      SafePutPixelSDL(screen, vertexPixels[i].x+1, vertexPixels[i].y+1, color);
      SafePutPixelSDL(screen, vertexPixels[i].x+1, vertexPixels[i].y, color);
      SafePutPixelSDL(screen, vertexPixels[i].x-1, vertexPixels[i].y-1, color);
      SafePutPixelSDL(screen, vertexPixels[i].x-1, vertexPixels[i].y+1, color);
      SafePutPixelSDL(screen, vertexPixels[i].x-1, vertexPixels[i].y, color);
    }
}
*/

void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result) {
    int N = result.size(); //we have to know the size in advance

    //Set initial values
    float currentX = a.x;
    float currentY = a.y;
    float currentZ = a.z;
    vec4 current3DPos = a.pos3d;

    //Set step values
    float x = (b.x - a.x) / float(max(N-1,1));
    float y = (b.y - a.y) / float(max(N-1,1));
    float z = (b.z - a.z) / float(max(N-1,1));
    vec4 step3d = (b.pos3d - a.pos3d) / float(max(N-1,1));

    //Calculate and store each step
    for(int i=0; i<N; i++){
        result[i].x = currentX;
        result[i].y = currentY;
        result[i].z = currentZ;
        result[i].pos3d = current3DPos;

        currentX = currentX+x;
        currentY = currentY+y;
        currentZ = currentZ+z;
        current3DPos = current3DPos+step3d;
    }
}

void SafePutPixelSDL(screen* screen, int x, int y, vec3 color) {
    if ((x >= 0) && (x < SCREEN_WIDTH) && (y >= 0) && (y < SCREEN_HEIGHT)){
        PutPixelSDL(screen, x, y, color);
    }
}

void TransformationMatrix(mat4& M, vec4 pos, vec4 rot) {

    mat4 toOrigin (1,0,0,-pos.x,
                    0,1,0,-pos.y,
                    0,0,1,-pos.z,
                    0,0,0,1);

    mat4 rotationX(1,0,0,0,
                    0,cos(rot.x),-sin(rot.x),0,
                    0,sin(rot.x),cos(rot.x),0,
                    0,0,0,1);

    mat4 rotationY(cos(rot.y),0,sin(rot.y),0,
                    0,1,0,0,
                    -sin(rot.y),0,cos(rot.y),0,
                    0,0,0,1);

    mat4 rotationZ(cos(rot.z),-sin(rot.z),0,0,
                    sin(rot.z),cos(rot.z),0,0,
                    0,0,1,0,
                    0,0,0,1);

    M = toOrigin*(rotationZ*rotationY*rotationX);
}

void VertexShader(const Vertex& v, Pixel& p){
    //input vertex v and set correct position and illumination for pixel p

    //Camera transform
    mat4 M;
    TransformationMatrix(M, cameraPos, cameraRot);
    vec4 localV = v.position*M; //for whatever reason they have to multiply this way around

    //Setting position
    p.z = 1.0f/glm::abs(glm::length(cameraPos-v.position));
    p.x = (focalLength * (localV.x/localV.z) + (SCREEN_WIDTH/2));
    p.y = (focalLength * (localV.y/localV.z) + (SCREEN_HEIGHT/2));
    p.pos3d = v.position;
}

vec3 getBand(float distance) {
    if(distance < band1.x) {
        return vec3(2);
    }
    if(distance < band2.x) {
            return vec3(0.5);
    }
    else {
        return vec3(0.25);
    }
}

void PixelShader(const Pixel& p, screen *screen, vec3 color, vec4 cNormal){
    //Draw the pixel p on the screen, if it is closest to the screen
    vec3 reflectance(1);
    int x = p.x;
    int y = p.y;
    if(p.z > depthBuffer[y][x]){
        depthBuffer[y][x] = p.z;

        float distance = glm::length((lightPos-p.pos3d));
        vec3 D = lightPower*max((float)0, glm::abs(glm::dot(cNormal,(lightPos-p.pos3d))));
        D = D*(float)(1/(4*glm::length(p.pos3d-lightPos)*M_PI));

        vec3 illumination = reflectance*(D + indirectLightPowerPerArea);
        vec3 band = getBand(distance);
        SafePutPixelSDL(screen, x, y, band*color);
    }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels){
    //given the vertecies of a triangle, produce two lists of pixels which
    //are at the left and rightmost positions on each row of the triangle

    //get maximum and minimum y values of the triangle
    int yMin = SCREEN_WIDTH, yMax = -1;
    for(int i = 0; i < 3; i++) {
        if(vertexPixels[i].y > yMax) {
            yMax = vertexPixels[i].y;
        }
        if(vertexPixels[i].y < yMin) {
            yMin = vertexPixels[i].y;
        }
    }
    //calculate number of rows needed to draw the triangle
    int ROWS = yMax - yMin + 1;
    leftPixels.resize(ROWS);
    rightPixels.resize(ROWS);

    //set default values for the vectors
    for(int i = 0; i < ROWS; i++) {
        //set left values as far right as possible and vice versa
        leftPixels[i].x = SCREEN_WIDTH;
        rightPixels[i].x = -1;

        //set y values to the height of the pixels that they represent
        leftPixels[i].y = yMin+i;
        rightPixels[i].y = yMin+i;
    }


    for (int q = 0; q < 3; q++){

        //hacky way of getting p and q to be every possible pair of vertices of the triangle
        int p = q + 1;
        if (q > 1){
          p = 0;
        }

        //get line between vertices p and q and ensure that it is long enough in y direction to have a pixel for every row
        int yDiff = glm::abs(vertexPixels[q].y - vertexPixels[p].y);
        vector<Pixel> line(yDiff+1);

        InterpolatePixel(vertexPixels[q], vertexPixels[p], line);

        //use line to set x values for leftPixels and rightPixels
        for (int i = 0; i < line.size(); i++) {
            for (int j = 0; j < leftPixels.size(); j++) {
                if (line[i].y == leftPixels[j].y){
                    if (line[i].x > rightPixels[j].x){
                        rightPixels[j].x = line[i].x;
                        rightPixels[j].z = line[i].z;
                        rightPixels[j].pos3d = line[i].pos3d;
                    }
                    if (line[i].x < leftPixels[j].x){
                        leftPixels[j].x = line[i].x;
                        leftPixels[j].z = line[i].z;
                        leftPixels[j].pos3d = line[i].pos3d;
                    }
                }
            }
        }

    }
}

void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, screen* screen, vec3 color, vec4 cNormal){
    //given left and right pixels of a triangle, call pixelshader for all pixels in the triangle
    for (int j = 0; j < leftPixels.size(); j++) {
        vector<Pixel> line(max(1,rightPixels[j].x - leftPixels[j].x +1)); // hacky shit
        InterpolatePixel(leftPixels[j], rightPixels[j], line);
        for(int i = 0; i < line.size(); i++) {
            if(line[i].y >= 0 && line[i].x >= 0 && line[i].x < SCREEN_WIDTH && line[i].y < SCREEN_HEIGHT) {
                PixelShader(line[i], screen, color, cNormal);
            }
        }
    }
}

void DrawPolygon( const vector<Vertex>& vertices, screen* screen, vec3 color, vec4 cNormal){
    //map vector of 4d points into vector of 2d points
    int vs = vertices.size();
    vector<Pixel> vertexPixels(vs);
    for(int i = 0; i < vs; i++) {
        VertexShader(vertices[i], vertexPixels[i]);
    }

    //fill leftPixels and rightPixels vectors
    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
    //draw the rows
    DrawRows(leftPixels, rightPixels, screen, color, cNormal);
}

void outlinePixel(screen* screen, int x, int y) {
    SafePutPixelSDL(screen, x-1, y-1, vec3(0,0,0));
    SafePutPixelSDL(screen, x-1, y, vec3(0,0,0));
    SafePutPixelSDL(screen, x-1, y+1, vec3(0,0,0));
    SafePutPixelSDL(screen, x, y-1, vec3(0,0,0));
    SafePutPixelSDL(screen, x, y, vec3(0,0,0));
    SafePutPixelSDL(screen, x, y+1, vec3(0,0,0));
    SafePutPixelSDL(screen, x+1, y-1, vec3(0,0,0));
    SafePutPixelSDL(screen, x+1, y, vec3(0,0,0));
    SafePutPixelSDL(screen, x+1, y+1, vec3(0,0,0));
}


/*Place your drawing here*/
void Draw(screen* screen, const vector <Triangle>& triangles){

    //reset depth buffer
    for( int y=0; y<SCREEN_HEIGHT; ++y ) {
        for( int x=0; x<SCREEN_WIDTH; ++x ) {
            depthBuffer[y][x] = 0;
        }
    }

    //reset screen
    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
    vec4 cNormal;
    //#pragma omp parallel for
    for( uint32_t i=0; i<triangles.size(); ++i ) {
        //set vertices array for a specific triangle
        vector<Vertex> vertices(3);
        vertices[0].position = triangles[i].v0;
        vertices[1].position = triangles[i].v1;
        vertices[2].position = triangles[i].v2;
        cNormal = triangles[i].normal;
        DrawPolygon(vertices, screen, triangles[i].color, cNormal);
    }
    //shader

    std::map<int,float> surround;
    int flag1 = 0;
    //#pragma omp parallel for
    for(int y = 1; y < SCREEN_HEIGHT-1; y++) {
        for(int x = 1; x < SCREEN_WIDTH-1; x++) {
            flag1 = 0;
            surround[1] = depthBuffer[y-1][x-1];
            surround[2] = depthBuffer[y-1][x];
            surround[3] = depthBuffer[y-1][x+1];
            surround[4] = depthBuffer[y][x-1];
            surround[5] = depthBuffer[y][x+1];
            surround[6] = depthBuffer[y+1][x-1];
            surround[7] = depthBuffer[y+1][x];
            surround[8] = depthBuffer[y+1][x+1];
            //ignore pixels outside cube
            for(int q = 1; q < 9; q++) {
                if(surround[q] == 0) {
                    flag1 = 1;
                }
            }
            if(flag1 != 1) {
                for(int z = 1; z < 9; z++) {
                    //if a surrounding pixel depth - main pixel depth is greater than certain distance
                    if(glm::abs(surround[z] - depthBuffer[y][x]) > 0.05f) {
                        outlinePixel(screen, x, y);
                        break;
                    }
                }
            }
        }
    }
}

/*Place updates of parameters here*/
void Update(){
    /* //Compute frame time
    static int t = SDL_GetTicks();
    int t2 = SDL_GetTicks();
    float dt = float(t2-t);
    t = t2;
    std::cout << "Render time: " << dt << " ms." << std::endl;
    //*/


    const uint8_t* keystate = SDL_GetKeyboardState( 0 );
    // lightPos.x = sin(theta) * 0.8;
    // lightPos.z = cos(theta) * 0.8;
    // if(keystate[SDL_SCANCODE_Z]){
    //     theta = theta + 0.1;
    // }
    // if(keystate[SDL_SCANCODE_X]){
    //     theta = theta - 0.1;
    // }


    float lookSpeed = 0.02;
    float moveSpeed = 0.02;

    //Collect button inputs
    vec4 lookVector(0, 0, 0, 1);
    vec4 moveVector(0, 0, 0, 1);

    if(keystate[SDL_SCANCODE_LEFT]){
        lookVector.y = -1;
    }
    else if(keystate[SDL_SCANCODE_RIGHT]){
        lookVector.y = 1;
    }
    if(keystate[SDL_SCANCODE_UP]){
        lookVector.x = 1;
    }
    else if(keystate[SDL_SCANCODE_DOWN]){
        lookVector.x = -1;
    }
    if(keystate[SDL_SCANCODE_W]){
        moveVector.z = -1;
    }
    else if(keystate[SDL_SCANCODE_S]){
        moveVector.z = 1;
    }
    if(keystate[SDL_SCANCODE_A]){
        moveVector.x = 1;
    }
    else if(keystate[SDL_SCANCODE_D]){
        moveVector.x = -1;
    }
    if(keystate[SDL_SCANCODE_SPACE]){
        moveVector.y = 1;
    }
    else if(keystate[SDL_SCANCODE_LCTRL]){
        moveVector.y = -1;
    }

    //*
    //Modify global variables
    mat4 rotate, rotation, translate;

    TransformationMatrix(rotate, lookVector*lookSpeed, vec4(0,0,0,0)); //translation of vector which stores camera rotation values
    TransformationMatrix(rotation, vec4(0,0,0,0), cameraRot*rotate); //rotation by new camera rotation values
    TransformationMatrix(translate, rotation*moveVector*moveSpeed, vec4(0,0,0,0)); //translation of vector which stores camera position values

    cameraRot = cameraRot*rotate;
    cameraPos = cameraPos*translate;
    //*/
}

int main(int argc, char* argv[]){
    screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
    vector<Triangle> triangles;
    LoadTestModel(triangles);

    time_t time0;
    time_t time1;
    omp_set_num_threads(4);
    time(&time0);
    int counter = 0;
    cout.precision(5);
    while(NoQuitMessageSDL()){
        Update();
        Draw(screen, triangles);
        SDL_Renderframe(screen);
        counter++;
        if(counter == 500) {
            time(&time1);

            float seconds = time1 - time0;
            int fps = 500/seconds;
            cout << fps << '\n';

            counter = 0;
            time(&time0);
        }
    }

    SDL_SaveImage(screen, "screenshot.bmp");
    KillSDL(screen);
    return 0;
}
