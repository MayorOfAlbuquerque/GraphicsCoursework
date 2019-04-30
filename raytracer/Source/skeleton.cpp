#include <iostream>
#include <omp.h>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

struct Intersection{
  vec4 position;
  float distance;
  int triangleIndex;
};

//lightSource with constants
vec4 lightPos;

#define SCREEN_WIDTH 512
#define SCREEN_HEIGHT 512
#define FULLSCREEN_MODE false

vec4 cameraPos(0, 0, -3.001, 1); //removing the 0.001 will cause a crash to occour
vec4 cameraRot(0, 0, 0, 1);
vec4 cameraDir(0, 0, 1, 0);

float theta = 0; //stores the rotation of the light around the room
bool blackWhiteMirror = false;

//float blur = 2;
float focus = 0;

float depBuf[SCREEN_HEIGHT][SCREEN_WIDTH];

vec3 colBuf[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 blurBuf[SCREEN_HEIGHT][SCREEN_WIDTH];

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void TransformationMatrix(mat4& M, vec4 pos, vec4 rot){
    mat4 toOrigin (1,0,0,-pos.x,
                    0,1,0,-pos.y,
                    0,0,1,-pos.z,
                    0,0,0,1);

    /* //don't seem to need this at the moment
    mat4 toCamera (1,0,0,pos.x,
                    0,1,0,pos.y,
                    0,0,1,pos.z,
                    0,0,0,1);
    //*/

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

void Update(){
    //Collect button inputs
    float lookSpeed = 0.02;
    float moveSpeed = 0.02;

    //float blurSpeed = 0.2;
    float focusSpeed = 0.1;

    vec4 lookVector(0, 0, 0, 1);
    vec4 moveVector(0, 0, 0, 1);
    const uint8_t* keystate = SDL_GetKeyboardState( 0 );

    if(keystate[SDL_SCANCODE_LSHIFT]){
      moveSpeed = 3*moveSpeed;
      lookSpeed = 3*lookSpeed;
    }

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

    /*
    if(keystate[SDL_SCANCODE_O]){
      if (blur <= (min(SCREEN_HEIGHT, SCREEN_WIDTH) - blurSpeed)){
        blur += blurSpeed;
      }
    }
    else if(keystate[SDL_SCANCODE_P]){
      if (blur > blurSpeed){
        blur -= blurSpeed;
      }
    }
    */

    if(keystate[SDL_SCANCODE_K]){
      focus += focusSpeed;
    }
    else if(keystate[SDL_SCANCODE_L]){
      focus -= focusSpeed;
    }

    //Modify global variables
    mat4 rotate, rotation, translate;

    TransformationMatrix(rotate, lookVector*lookSpeed, vec4(0)); //translation of vector which stores camera rotation values
    TransformationMatrix(rotation, vec4(0), cameraRot*rotate); //rotation by new camera rotation values
    TransformationMatrix(translate, rotation*moveVector*moveSpeed, vec4(0)); //translation of vector which stores camera position values

    cameraRot = cameraRot*rotate;
    cameraPos = cameraPos*translate;
}

bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection){
    bool intersection = false;

    for (uint32_t i=0; i < triangles.size(); i++){
    Triangle triangle = triangles[i];

    vec4 v0 = triangle.v0;
    vec4 v1 = triangle.v1;
    vec4 v2 = triangle.v2;

    vec3 e1 = vec3(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);
    vec3 e2 = vec3(v2.x-v0.x, v2.y-v0.y, v2.z-v0.z);
    vec3 b = vec3(start.x-v0.x, start.y-v0.y, start.z-v0.z);
    mat3 A(-vec3(dir), e1, e2);
    vec3 x = glm::inverse(A) * b;

    float t = x.x;
    float u = x.y;
    float v = x.z;

        if ((t >= 0) && (u >= 0) && (v >= 0) && (u + v <= 1)){ //use less than or equal to here, instead of less than
            intersection = true;

            if (t <= closestIntersection.distance){
            closestIntersection.position = (start + (dir * t));
            closestIntersection.distance = t;
            closestIntersection.triangleIndex = i;
            }
        }
    }

    return intersection;
}

vec3 DirectLight(const Intersection& i, vec4 lightPos, vec3 lightColor, const vector<Triangle>& triangles) {
    vec3 bwColour;

    vec4 dirToLight = glm::normalize(lightPos - i.position);
    float disToLight = glm::length(lightPos - i.position);

    Intersection closestIntersection = {i.position, std::numeric_limits<float>::max(), -1};
    ClosestIntersection(i.position+0.01f*dirToLight, dirToLight, triangles, closestIntersection);

    if ((disToLight <= closestIntersection.distance)){ //if not in shadow, add direct light
        bwColour = (float) (1.0f/(4.0f * M_PI * disToLight)) * lightColor;
    }

    bwColour += vec3(0.2f); //add ambient light

    return bwColour;
}

void godBlur(float a, float minDepth, float maxDepth){
    int r = (int) a;

    minDepth += focus;
    maxDepth += focus;

    if (r < 1){ //ensure that r is within the correct bounds
        r = 1;
    }
    else if (r > min(SCREEN_WIDTH, SCREEN_HEIGHT)){
        r = min(SCREEN_WIDTH, SCREEN_HEIGHT);
    }
    if (r%2 != 1){ //ensure that r is odd
        r -= 1;
    }

    int neigh = (r-1)/2; //number of neighbours on either side of the current pixel;

    //perform horizontal blur
    for (int i = 0; i < SCREEN_HEIGHT; i++){ //loop over rows

        vec3 acc = vec3(0); //accumulator for row
        int li = 0; //index for left part of row
        int ri = neigh; //index for right part of row
        int ti = neigh + 1; //total number of cells in accumulator

        for (int j = 0; j <= ri; j++){ //initialize accumulator
            acc = acc + colBuf[i][j];
        }
        //*
        for (int j = 0; j < neigh; j++){ //blur cells with not enough neighbours on the left
            blurBuf[i][j] = acc/(float) ti;
            ri += 1;
            ti += 1;
            acc = acc + colBuf[i][ri];
        }
        //*/
        //*
        for (int j = neigh; j < SCREEN_WIDTH - neigh; j++){ //blur cells with enough neighbours on both sides
            blurBuf[i][j] = acc/(float) ti;
            acc = acc - colBuf[i][li];
            li += 1;
            ri += 1;
            acc = acc + colBuf[i][ri];
        }
        //*/
        //*
        for (int j = SCREEN_WIDTH - neigh; j < SCREEN_WIDTH; j++){ //blur cells with not enough neighbours on the right
            blurBuf[i][j] = acc/(float) ti;
            acc = acc - colBuf[i][li];
            li += 1;
            ti -= 1;
        }
        //*/
    }

    //perform vertical blur
    for (int i = 0; i < SCREEN_WIDTH; i++){ //loop over columns

        /*
        for (int j = 0; j < SCREEN_HEIGHT; j++){
        colBuf[j][i] = blurBuf[j][i];
        }
        //*/

        vec3 acc = vec3(0); //accumulator for column
        int li = 0; //index for top part of column
        int ri = neigh; //index for bottom part of column
        int ti = neigh + 1; //total number of cells in accumulator

        for (int j = 0; j <= ri; j++){ //initialize accumulator
            acc = acc + blurBuf[j][i];
        }
        //*
        for (int j = 0; j < neigh; j++){ //blur cells with not enough neighbours on the left
            if ((depBuf[j][i] >= minDepth) && (depBuf[j][i] <= maxDepth)){
                colBuf[j][i] = acc/(float) ti;
            }
        ri += 1;
        ti += 1;
        acc = acc + blurBuf[ri][i];
        }
        //*/
        //*
        for (int j = neigh; j < SCREEN_HEIGHT - neigh; j++){ //blur cells with enough neighbours on both sides
            if ((depBuf[j][i] >= minDepth) && (depBuf[j][i] <= maxDepth)){
                colBuf[j][i] = acc/(float) ti;
            }
            acc = acc - blurBuf[li][i];
            li += 1;
            ri += 1;
            acc = acc + blurBuf[ri][i];
        }
        //*/
        //* problem is happening here somehow
        for (int j = SCREEN_HEIGHT - neigh; j < SCREEN_HEIGHT; j++){ //blur cells with not enough neighbours on the right
            if ((depBuf[j][i] >= minDepth) && (depBuf[j][i] <= maxDepth)){
                colBuf[j][i] = acc/(float) ti;
            }
            acc = acc - blurBuf[li][i];
            li += 1;
            ti -= 1;
        }
        //*/
    }
}

// vec3 SoftLight(const Intersection& i, vec4 lightPos, vec3 lightColor, const vector<Triangle>& triangles) {
//     vec3 bwColour = vec3(0);
//     for(int x = 0; x < 20; x++) {
//         if(bwColour.x > 1 || bwColour.y > 1 || bwColour.z > 1) break;
//         for(int y = 0; y < 20; y++) {
//             if(bwColour.x > 1 || bwColour.y > 1 || bwColour.z > 1) break;
//             lightPos.z = lightPos.z+0.1;
//             bwColour += DirectLight(i, lightPos, lightColor, triangles);
//         }
//         lightPos.x = lightPos.x+0.1;
//     }
//     return bwColour;
// }

void DrawMirroredWall(int x, int y, vec4 lightPos, vec3 lightColor, vec4 dir, Triangle tri, Intersection closestIntersection, const vector<Triangle>& triangles, screen* screen) {
    //get angle of line
    //get angle of line after reflection
    vec4 reflectedLine = glm::reflect(dir, tri.normal);
    //cout << reflectedLine.z << "\n";
    vec4 mirrorIntersection = closestIntersection.position;
    vec4 norm = tri.normal;
    vec3 mirrTint = tri.color;
    mirrorIntersection = mirrorIntersection + norm * vec4(0.001f);
    //find object hit next
    closestIntersection = {mirrorIntersection, std::numeric_limits<float>::max(), -1};
    //colour this pixel that colour - add shadows later
    if (ClosestIntersection(mirrorIntersection, reflectedLine, triangles, closestIntersection)){
       
        //get colour of reflected object
        vec3 colour = triangles[closestIntersection.triangleIndex].color;
        //get average of colour values
        vec3 bwColour = DirectLight(closestIntersection, lightPos, lightColor, triangles);
        if(blackWhiteMirror) {
            float avColourVal = (colour.x + colour.y + colour.z)/3;
            colBuf[y][x] = bwColour*avColourVal;
            //PutPixelSDL(screen, x, y, avColourVal*bwColour);
        }
        else {
            colBuf[y][x] = bwColour*colour*mirrTint;
            //PutPixelSDL(screen, x, y, colour*bwColour);
        }
       
    }
    else {
        //put black pixel if it points to empty space
        colBuf[y][x] = mirrTint*vec3(0.05);
        //PutPixelSDL(screen, x, y, vec3(0));
    }
}



void Draw(screen* screen, const vector<Triangle>& triangles){

    //Clear the buffer
    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

    //since we're using xf and yf, focal length needs to be equal to width of screen which is the magnitude of the interval [-1:1] which is 2
    float focalLength = 1;

    vec3 lightColor = 14.f * vec3(1);

    //lightSource rotates around the room
    // lightPos.x = sin(theta) * 0.8;
    // lightPos.z = cos(theta) * 0.8;
    // theta = theta + 0.05;

    //Loop over all pixels in the image and calculate a ray for each one.
    #pragma omp parallel for
    for (int y = 0; y < SCREEN_HEIGHT; y++){ //don't use unsigned ints here!!!
        for (int x = 0; x < SCREEN_WIDTH; x++){
            //get direction to pixel from middle of camera view
            float xf = (float) x / (SCREEN_WIDTH-1) - 0.5; //goes from interval [0:SCREEN_WIDTH-1] to [-1/2:1/2]
            float yf = (float) y / (SCREEN_HEIGHT-1) - 0.5; //goes from interval [0:SCREEN_HEIGHT-1] to [-1/2:1/2]
            vec4 dir(xf, yf, focalLength, 0);

            //modify direction to pixel based on camera rotation
            mat4 M;
            TransformationMatrix(M, vec4(0), cameraRot);
            dir = M * dir;

            

            Intersection closestIntersection = {cameraPos, std::numeric_limits<float>::max(), -1};
            if (ClosestIntersection(cameraPos, dir, triangles, closestIntersection)){
                Triangle tri = triangles[closestIntersection.triangleIndex];
                if(tri.mirror) {
                DrawMirroredWall(x, y, lightPos, lightColor, dir, tri, closestIntersection, triangles, screen);
                }
                else {
                    depBuf[y][x] = closestIntersection.distance;
                    vec3 bwColour = DirectLight(closestIntersection, lightPos, lightColor, triangles);
                    vec3 colour = triangles[closestIntersection.triangleIndex].color;
                    colBuf[y][x] = bwColour*colour;
                }
            }
            else{
                colBuf[y][x] = vec3(0);
                depBuf[y][x] = std::numeric_limits<float>::max();
            }
        }
    }

    godBlur(11, -100.0, 1.0); //apply blur
    godBlur(9, 1.0, 2.0); //apply blur
    godBlur(7, 2.0, 3.0); //apply blur
    godBlur(5, 3.0, 4.0); //apply blur
    godBlur(3, 4.0, 5.0); //apply blur
    //godBlur(9, 5.0, 6.0); //apply blur
    godBlur(3, 6.0, 7.0); //apply blur
    godBlur(5, 7.0, 8.0); //apply blur
    godBlur(7, 8.0, 9.0); //apply blur
    godBlur(9, 9.0, 10.0); //apply blur
    godBlur(11, 10.0, 100.0); //apply blur

    for (int y = 0; y < SCREEN_HEIGHT; y++){ //don't use unsigned ints here!!!
        for (int x = 0; x < SCREEN_WIDTH; x++){
          PutPixelSDL(screen, x, y, colBuf[y][x]);
        }
    }
}

int main(int argc, char* argv[]){
    //Initialize the screen
    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

    //Load the test model into an empty vector of triangles
    vector<Triangle> triangles;
    LoadTestModel(triangles);

    lightPos.y = -0.5;
    lightPos.w = 1.0;
    //Enter the rendering loop
    while(NoQuitMessageSDL()){
        cout << "Focus: " << focus << "\n";
        Update();
        Draw(screen, triangles);
        SDL_Renderframe(screen);
    }

    //Finalise
    SDL_SaveImage( screen, "screenshot.bmp" );
    KillSDL(screen);
    return 0;
}
