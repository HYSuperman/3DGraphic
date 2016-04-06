////////////////////////////////////////////////////////////////////////
//
//   Harvard University
//   CS175 : Computer Graphics
//   Professor Steven Gortler
//
////////////////////////////////////////////////////////////////////////
//	These skeleton codes are later altered by Ming Jin,
//	for "CS6533: Interactive Computer Graphics", 
//	taught by Prof. Andy Nealen at NYU
////////////////////////////////////////////////////////////////////////

#include <list>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <algorithm>

// #if __GNUG__
// #   include <tr1/memory>
// #endif

#include <GL/glew.h>
#ifdef __MAC__
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

#include "cvec.h"
#include "matrix4.h"
#include "ppm.h"
#include "glsupport.h"

// hw3 include ------------
#include "arcball.h"
#include "geometrymaker.h"
#include "rigtform.h"
// ------------------------

// hw4 include ------------
#include "asstcommon.h"
#include "scenegraph.h"
#include "drawer.h"
#include "picker.h"
// ------------------------

// hw5 include ------------
#include "sgutils.h"
// ------------------------

using namespace std;      // for string, vector, iostream, and other standard C++ stuff
// using namespace tr1; // for shared_ptr

// G L O B A L S ///////////////////////////////////////////////////

// --------- IMPORTANT --------------------------------------------------------
// Before you start working on this assignment, set the following variable
// properly to indicate whether you want to use OpenGL 2.x with GLSL 1.1 or
// OpenGL 3.x+ with GLSL 1.3.
//
// Set g_Gl2Compatible = true to use GLSL 1.0 and g_Gl2Compatible = false to
// use GLSL 1.3. Make sure that your machine supports the version of GLSL you
// are using. In particular, on Mac OS X currently there is no way of using
// OpenGL 3.x with GLSL 1.3 when GLUT is used.
//
// If g_Gl2Compatible=true, shaders with -gl2 suffix will be loaded.
// If g_Gl2Compatible=false, shaders with -gl3 suffix will be loaded.
// To complete the assignment you only need to edit the shader files that get
// loaded
// ----------------------------------------------------------------------------

const bool g_Gl2Compatible = true; // hw4

static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static const float g_frustNear = -0.1;    // near plane
static const float g_frustFar = -50.0;    // far plane
static const float g_groundY = -2.0;      // y coordinate of the ground
static const float g_groundSize = 10.0;   // half the ground length

static int g_windowWidth = 512;
static int g_windowHeight = 512;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;


// hw2 globals -------------------------------------------------
// number of cubes to draw
static const int CUBE_NUMBER = 2;
// number of robots to draw
static const int ROBOT_NUMBER = 2;
// default, sky camera frame
// sky camera frame: 0; robot one: 1; robot two: 2;
static int CHOOSEN_FRAME = 0;
// default, world sky camera
static bool WORLD_SKY = true;
// just names
static string MODE[] = {"sky", "robot 1", "robot 2"};
static string SKY_MODE[] = {"world-sky", "sky-sky"};
//  -------------------------------------------------------------


// hw3 globals --------------------------------------------------
// default radius on the screen
static double g_arcballScreenRadius = 0.25 * min(g_windowWidth, g_windowHeight);
// scale between obj frame and screen. ?
static double g_arcballScale = getScreenToEyeScale(-4.0, g_frustFovY, g_windowHeight);

// consts about the sphere
static const float sphere_radius = 1.0;
static const int sphere_slice = 30;
static const int sphere_stack = 30;
static bool DRAW_SPHERE = true;

// --------------------------------------------------------------



// picking shaders

static bool PICKING = false;
static const int PICKING_SHADER = 2; // index of the picking shader is g_shaerFiles
static const int g_numShaders = 3;
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader"}
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
  {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/pick-gl2.fshader"}
};
static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states


// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// A vertex with floating point position and normal
struct VertexPN {
  Cvec3f p, n;

  VertexPN() {}
  VertexPN(float x, float y, float z,
           float nx, float ny, float nz)
    : p(x,y,z), n(nx, ny, nz)
  {}

  // Define copy constructor and assignment operator from GenericVertex so we can
  // use make* functions from geometrymaker.h
  VertexPN(const GenericVertex& v) {
    *this = v;
  }

  VertexPN& operator = (const GenericVertex& v) {
    p = v.pos;
    n = v.normal;
    return *this;
  }
};

struct Geometry {
  GlBufferObject vbo, ibo;
  int vboLen, iboLen;

  Geometry(VertexPN *vtx, unsigned short *idx, int vboLen, int iboLen) {
    this->vboLen = vboLen;
    this->iboLen = iboLen;

    // Now create the VBO and IBO
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPN) * vboLen, vtx, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
  }

  void draw(const ShaderState& curSS) {
    // Enable the attributes used by our shader
    safe_glEnableVertexAttribArray(curSS.h_aPosition);
    safe_glEnableVertexAttribArray(curSS.h_aNormal);

    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, p));
    safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, n));

    // bind ibo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

    // draw!
    glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);

    // Disable the attributes used by our shader
    safe_glDisableVertexAttribArray(curSS.h_aPosition);
    safe_glDisableVertexAttribArray(curSS.h_aNormal);
  }
};

typedef SgGeometryShapeNode<Geometry> MyShapeNode; // hw4


// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube;
// hw3 a shared sphere
static shared_ptr<Geometry> g_sphere;

// --------- Scene

static shared_ptr<SgRootNode> g_world;
static shared_ptr<SgRbtNode> g_skyNode, g_groundNode, g_robot1Node, g_robot2Node;
static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do picking
static shared_ptr<SgRbtNode> g_currentView;

static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
static RigTForm g_skyRbt(Cvec3(0.0, 0.25, 4.0)); // hw3 RigTForm g_skyRbt

// hw3 global
// sphere center
static Cvec3 sphere_center_eyeCord; 

// function signature for reset() use
static void initScene();


static list<vector<RigTForm> > keyFrames;
static list<vector<RigTForm> >::iterator currentKeyFrame = keyFrames.end(); 
static int currentIndex = -1;
static vector<shared_ptr<SgRbtNode> > nodeData;

static int g_msBetweenKeyFrames = 2000; // 2 seconds between keyframes
static int g_animateFramesPerSecond = 60; // frames to render per second during animation playback
static bool PLAYING = false;


///////////////// END OF GLOBALS ///////////////////


static void initGround() {
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
  VertexPN vtx[4] = {
    VertexPN(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
    VertexPN(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
  };
  unsigned short idx[] = {0, 1, 2, 0, 2, 3};
  g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6));
}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

// hw3 init a sphere
static void initSphere(){
  int ibLen, vbLen;
  getSphereVbIbLen(sphere_slice, sphere_stack, vbLen, ibLen);

  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);
  
  makeSphere(sphere_radius, sphere_slice, sphere_stack, vtx.begin(), idx.begin());
  g_sphere.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
  GLfloat glmatrix[16];
  projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
  safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
}

// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
static void updateFrustFovY() {
  if (g_windowWidth >= g_windowHeight)
    g_frustFovY = g_frustMinFov;
  else {
    const double RAD_PER_DEG = 0.5 * CS175_PI/180;
    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
           g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}

static void drawStuff(const ShaderState& curSS, bool picking) {
  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  // Choose an eye camera
  // changed to shared_ptr in hw4
  switch(CHOOSEN_FRAME){
    case 0:
      g_currentView = g_skyNode;
      break;
    case 1:
      g_currentView = g_robot1Node;
      break;
    case 2:
      g_currentView = g_robot2Node;
      break;
  }

  const RigTForm eyeRbt = getPathAccumRbt(g_world, g_currentView);
  const RigTForm invEyeRbt = inv(eyeRbt);

  const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
  safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
  safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

  if (!picking) {
    Drawer drawer(invEyeRbt, curSS);
    g_world->accept(drawer);

    // draw a sphere?
    if(!DRAW_SPHERE){return;}    
    if(!WORLD_SKY && !g_currentPickedRbtNode){return;}
    
    RigTForm sphere_center = getPathAccumRbt(g_world, g_currentPickedRbtNode);
    RigTForm eyeCord = invEyeRbt * sphere_center; 
    sphere_center_eyeCord = eyeCord.getTranslation();

    if (!(g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton))){ 
      g_arcballScale = getScreenToEyeScale(sphere_center_eyeCord[2], g_frustFovY, g_windowHeight);
    }
    double arcball_radius = g_arcballScale * g_arcballScreenRadius;

    // draw sphere
    // ----------
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    Matrix4 MVM = rigTFormToMatrix(eyeCord) * Matrix4::makeScale(Cvec3(arcball_radius, arcball_radius, arcball_radius));
    Matrix4 NMVM = normalMatrix(MVM);
    sendModelViewNormalMatrix(curSS, MVM, NMVM);
    safe_glUniform3f(curSS.h_uColor, 0, 0, 0);
    g_sphere->draw(curSS);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  }
  else{
    Picker picker(invEyeRbt, curSS);
    g_world->accept(picker);
    glFlush();
    g_currentPickedRbtNode = picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
    if (g_currentPickedRbtNode == g_groundNode)
      g_currentPickedRbtNode = shared_ptr<SgRbtNode>(); 
  }  

}

static void pick() {
  // We need to set the clear color to black, for pick rendering.
  // so let's save the clear color
  GLdouble clearColor[4];
  glGetDoublev(GL_COLOR_CLEAR_VALUE, clearColor);

  glClearColor(0, 0, 0, 0);

  // using PICKING_SHADER as the shader
  glUseProgram(g_shaderStates[PICKING_SHADER]->program);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawStuff(*g_shaderStates[PICKING_SHADER], true);

  // Uncomment below and comment out the glutPostRedisplay in mouse(...) call back
  // to see result of the pick rendering pass
  // glutSwapBuffers();

  //Now set back the clear color
  glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);

  checkGlErrors();
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  // IF in PICKING mode, use PICKING shader
  if(PICKING){
    pick();
    PICKING = false;
    cout << "PICKING: OFF" << endl;

  }
  else{
    drawStuff(*g_shaderStates[g_activeShader], false);
    glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)
  }
  checkGlErrors();
}

// callback when we change the size of window
static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  // hw3
  g_arcballScreenRadius = 0.25 * min(g_windowWidth, g_windowHeight);
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  updateFrustFovY();
  glutPostRedisplay();
}

// compute Quateinion matrix
// used below in motion()
static Quat computeSphereQuat(const int x, const int y, const double dx, const double dy){

  // hw3 computer on-screen sphere center
  Cvec2 sphere_center_screen = getScreenSpaceCoord(sphere_center_eyeCord, makeProjectionMatrix(), g_frustNear, g_frustFovY,
                                                   g_windowWidth, g_windowHeight);

  const double Dx0 = g_mouseClickX - sphere_center_screen[0], Dy0 = g_mouseClickY - sphere_center_screen[1];
  const double Dx1 = x - sphere_center_screen[0], Dy1 = g_windowHeight - y - 1 - sphere_center_screen[1];
  const double xyRadius12 = Dx0*Dx0 + Dy0*Dy0;
  const double xyRadius22 = Dx1*Dx1 + Dy1*Dy1;
  double z0 = 0.0, z1 = 0.0;

  if(g_arcballScreenRadius*g_arcballScreenRadius >= xyRadius12){
    z0 = sqrt(g_arcballScreenRadius * g_arcballScreenRadius - xyRadius12);
  }

  if(g_arcballScreenRadius*g_arcballScreenRadius >= xyRadius22){
    z1 = sqrt(g_arcballScreenRadius * g_arcballScreenRadius - xyRadius22);
  }

  Cvec3 v0(Dx0, Dy0, z0);
  Cvec3 v1(Dx1, Dy1, z1);

  // if there is no move, just return identity quat
  if(v0[0] == v1[0] && v0[1] == v1[1] && v0[2] == v1[2]){return Quat();}

  // cross product
  Cvec3 k = normalize(cross(v0, v1));
  // because it is arcball, we spin twice the angle!
  double angle_cos = dot(v0, v1) / (norm(v0)*norm(v1));
  double angle_sin = sqrt(1 - angle_cos*angle_cos);

  return Quat(angle_cos, k*angle_sin);
}

// compute AuxiliaryMatrix
// used below in motion()
static RigTForm computeAuxiliaryMatrix(){
  RigTForm objectFrame;
  RigTForm eyeFrame;
  RigTForm baseFrame;

  // objectFrame and eyeFrame
  objectFrame = getPathAccumRbt(g_world, g_currentPickedRbtNode);
  eyeFrame = getPathAccumRbt(g_world, g_currentView);

  // compute baseFrame
  baseFrame = getPathAccumRbt(g_world, g_currentPickedRbtNode, 1);

  // As = S^-1 * A
  return inv(baseFrame) * makeMixedFrame(objectFrame, eyeFrame);
}

// compute rbt matrix
// used below in motion()
static RigTForm computeRBT(Quat theQuat, const double dx, const double dy){

  RigTForm Q;
  
  if (g_mouseLClickButton && !g_mouseRClickButton) { // left button down?
    Q = RigTForm(theQuat);
  }
  else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
    if(DRAW_SPHERE)
      Q = RigTForm(Cvec3(dx, dy, 0) * g_arcballScale);//hw3
    else
      Q = RigTForm(Cvec3(dx, dy, 0) * 0.01); 
    if(!WORLD_SKY && !CHOOSEN_FRAME){Q = inv(Q);}
  }
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
    if(DRAW_SPHERE)
      Q = RigTForm(Cvec3(0, 0, -dy) * g_arcballScale);//hw3
    else
      Q = RigTForm(Cvec3(0, 0, -dy) * 0.01);
    if(!WORLD_SKY && !CHOOSEN_FRAME){Q = inv(Q);}
  }

  return Q;
}

// callback when we click mouse and move mouse
static void motion(const int x, const int y) {

  // movement on the screen
  const double dx = x - g_mouseClickX;
  const double dy = g_windowHeight - y - 1 - g_mouseClickY;

  // compute quaternion matrix
  Quat theQuat = computeSphereQuat(x, y, dx, dy);
  // compute RBT matrix
  RigTForm Q = computeRBT(theQuat, dx, dy);
  // compute auxiliary matrix
  RigTForm A = computeAuxiliaryMatrix();

  // change the nodes
  if (g_mouseClickDown) {

    if(!g_currentPickedRbtNode){
      if(WORLD_SKY) { g_skyNode->setRbt(doQtoOwrtA(inv(Q), g_skyNode->getRbt(), A));}
      else { g_skyNode->setRbt(doQtoOwrtA(inv(Q), g_skyNode->getRbt(), getPathAccumRbt(g_world, g_skyNode)));}
    }
    else{
      g_currentPickedRbtNode->setRbt(doQtoOwrtA(Q, g_currentPickedRbtNode->getRbt(), A));
    }
    glutPostRedisplay(); // we always redraw if we changed the scene
  }

  // update mouse position
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}

static void reset(){
  //reset g_skyRbt and g_objectRbt to their default values
  g_skyNode->setRbt(RigTForm(Cvec3(0.0, 0.25, 4.0)));
  // robots (hw4 and after)
  g_robot1Node->setRbt(RigTForm(Cvec3(-2, 1, 0)));
  g_robot2Node->setRbt(RigTForm(Cvec3(2, 1, 0)));
  g_currentView = shared_ptr<SgRbtNode>(); 
  g_currentPickedRbtNode = shared_ptr<SgRbtNode>();

  keyFrames = list<vector<RigTForm> >();
  currentKeyFrame = keyFrames.end(); 
  nodeData = vector<shared_ptr<SgRbtNode> >();

  //reset the views and manipulation mode to default
  CHOOSEN_FRAME = 0;
  //reset sky camera mode to use the "world-sky" frame
  WORLD_SKY = true;
  // draw sphere
  DRAW_SPHERE = true;

  PICKING = false;

	cout << "reset objects and modes to defaults" << endl;
  cout << "reset to draw sphere mode" << endl;

  initScene();
}

// callback when we click mouse
static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;

  // hw3
  // glutPostRedisplay();
  glutPostRedisplay();

}

static void frameToSgRbtNodes(vector<RigTForm> frame){

  vector<RigTForm>::iterator i;
  vector<shared_ptr<SgRbtNode> >::iterator j;

  for(i = frame.begin(), j = nodeData.begin(); i != frame.end() && j != nodeData.end(); i++,j++){
    (*j)->setRbt(*i);
  }

}

static vector<RigTForm> sgRbtNodesToFrame(bool newFrame){

  vector<RigTForm>::iterator i;
  vector<shared_ptr<SgRbtNode> >::iterator j;

  if(newFrame){
    vector<RigTForm> frame = vector<RigTForm>(nodeData.size(), RigTForm());
    for(i = frame.begin(), j = nodeData.begin(); i != frame.end() && j != nodeData.end(); i++,j++){
      (*i) = (*j)->getRbt();
    }
    return frame;
  }
  else{
    for(i = currentKeyFrame->begin(), j = nodeData.begin(); i != currentKeyFrame->end() && j != nodeData.end(); i++,j++){
      (*i) = (*j)->getRbt();
    }
    return (*currentKeyFrame);
  }

}


bool interpolateAndDisplay(float t) {

  int keyFrameSize = keyFrames.size();

  int before = (int)floor(t);
  int after = before + 1;
  double alpha = t - floor(t);

  if(before >= keyFrameSize - 3){
    list<vector<RigTForm> >::iterator t = keyFrames.end();
    t--;t--;
    currentKeyFrame = t;
    frameToSgRbtNodes(*currentKeyFrame);
    glutPostRedisplay();
    return true;
  }
  else{
    int index = before + 1;
    list<vector<RigTForm> >::iterator t = keyFrames.begin();
    while(index--){
      t++;
    }

    vector<RigTForm> showFrame_before = (*t);
    vector<RigTForm> showFrame_after = (*next(t));
    vector<RigTForm> showFrame;

    vector<RigTForm>::iterator i, j;
    for(i = showFrame_before.begin(), j = showFrame_after.begin();i != showFrame_before.end() && j != showFrame_after.end();i++, j++){
      showFrame.push_back(lerp((*i), (*j), alpha));
    }
    frameToSgRbtNodes(showFrame);
    glutPostRedisplay();
    return false;
  }

}

static void animateTimerCallback(int ms) {
  if(!PLAYING){
    cout << "Stopping Animation Playback...\n";
    return;
  }
  float t = (float)ms/(float)g_msBetweenKeyFrames;
  bool endReached = interpolateAndDisplay(t);
  if (!endReached){
    glutTimerFunc(1000/g_animateFramesPerSecond, animateTimerCallback, ms + 1000/g_animateFramesPerSecond);
  }
  else{
    PLAYING = false;
    cout << "You have reached the end of the animation, now showing the last frame.\n";
    return;
  } 
}

// callback when we press a key
static void keyboard(const unsigned char key, const int x, const int y) {

  switch (key) {
  case 27:
    exit(0);                                  // ESC

  // Spacebar
  case 32:
    if(currentKeyFrame != keyFrames.end()){
      cout << "Restoring to Current Frame: keyFrames[" << currentIndex << "].\n"; 
      frameToSgRbtNodes(*currentKeyFrame);
    }
    break;

  // -'a' control whether draw the sphere
  case 'a':
    DRAW_SPHERE = !DRAW_SPHERE;

    cout << "Now in ";
    if(DRAW_SPHERE)cout << "SPHERE";
    else cout << "NO-SPHERE";  
    cout << " mode.\n";
    break;


  case 'd':
    if(currentKeyFrame != keyFrames.end()){
      cout << "Deleting Current Frame: keyFrames[" << currentIndex << "].\n";
      // delete current key frame
      list<vector<RigTForm> >::iterator tempFrame = currentKeyFrame;
      tempFrame++;
      keyFrames.erase(currentKeyFrame);

      if(!keyFrames.empty()){
        if(tempFrame != keyFrames.begin()){
          currentKeyFrame = --tempFrame;
          currentIndex--;
        }
        else{
          currentKeyFrame = tempFrame;
        }
        cout << "Now Current Frame: keyFrames[" << currentIndex << "].\n";
        frameToSgRbtNodes((*currentKeyFrame));
      }
      else{
        cout << "Now it is an empty list.\n";
        currentKeyFrame = keyFrames.end();
        currentIndex = -1;
      }

    }
    break;
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "p\t\tChoose object to edit\n"
    << "v\t\tChoose view\n"
    << "drag left mouse to rotate\n" << endl;
    break;

  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;

  case 'f':
    g_activeShader ^= 1;
    break;

  // update a frame
  case 'u':
    // if not empty
    if(currentKeyFrame != keyFrames.end()){
      cout << "Updating Current Frame: keyFrames[" << currentIndex << "].\n";
      sgRbtNodesToFrame(false);
      break;
    }

  // add a new frame
  case 'n':
  {

    vector<RigTForm> newFrame = sgRbtNodesToFrame(true);
    // if not empty
    if(currentKeyFrame != keyFrames.end()){
      keyFrames.insert(next(currentKeyFrame), newFrame);
      currentKeyFrame++;
    }
    // if empty
    else{
      keyFrames.push_back(newFrame);
      currentKeyFrame = keyFrames.begin();
    }

    cout << "Added a new frame at keyFrames[" << ++currentIndex << "].\n";

    break;
  }

  case '>':
    if(currentKeyFrame != prev(keyFrames.end()) && currentKeyFrame != keyFrames.end()){
      currentKeyFrame++;
      currentIndex++;
      cout << "Forward to: keyFrames[" << currentIndex << "].\n";
      frameToSgRbtNodes((*currentKeyFrame));
    }
    break;

  case '<':
    if(currentKeyFrame != keyFrames.begin() && currentKeyFrame != keyFrames.end()){
      currentKeyFrame--;
      currentIndex--;
      cout << "Backward to: keyFrames[" << currentIndex << "].\n";
      frameToSgRbtNodes(*currentKeyFrame);
    }
    break;

  case 'i':
  {
    ifstream input;
    input.open("Frames.txt", ifstream::in);
    int frameNumber, frameSize;
    input >> frameNumber >> frameSize;
    keyFrames.clear();

    cout << "Reading keyFrameList from file Frame.txt.\n";
    for(int i = 0;i < frameNumber;i++){
      vector<RigTForm> tempFrame;
      for(int j = 0;j < frameSize;j++){
        Cvec3 t;
        Quat q;
        for(int k = 0;k < 7;k++){
          double temp;
          input >> temp;
          if(k < 3){t[k] = temp;}
          else{q[k - 3] = temp;}
        }
        tempFrame.push_back(RigTForm(t, q));
      }
      keyFrames.push_back(tempFrame);
    }

    currentKeyFrame = keyFrames.begin();
    if(keyFrames.empty()){
      currentIndex = -1;
    }
    else{
      currentIndex = 0;
    }
    frameToSgRbtNodes((*currentKeyFrame));
    input.close();

    break;
  }

  case 'w':
  {
    ofstream output;
    output.open("Frames.txt", ofstream::out);
    int frameNumber = keyFrames.size();
    int frameSize;

    if(!frameNumber){
      frameSize = 0;
    }
    else{
      frameSize = keyFrames.front().size();
    }

    cout << "Writring keyFrameList to file Frame.txt.\n";
    output << frameNumber << " " << frameSize << endl;
    for(list<vector<RigTForm> >::iterator i = keyFrames.begin();i != keyFrames.end();i++){
      for(vector<RigTForm>::iterator j = (*i).begin();j != (*i).end();j++){
        for(int i = 0;i < 3;i++){
          output << (*j).getTranslation()[i] << " ";
        }
        for(int i = 0;i < 4;i++){
          if(i < 3){output << (*j).getRotation()[i] << " ";}
          else{output << (*j).getRotation()[i] << endl;}
        }
      }
    }

    break;
  }

  case 'y':
    if(keyFrames.size() < 4){
      cout << "Warning: not enough frames, must be at lease 4.\n";
      break;
    }
    else{
      PLAYING = !PLAYING;
      if(PLAYING){
        cout << "Staring Animation Playback...\n";
        glutTimerFunc(0, animateTimerCallback, 0);
      }
      else{
        list<vector<RigTForm> >::iterator t = keyFrames.end();
        t--;t--;
        currentKeyFrame = t;
        currentIndex = keyFrames.size() - 2;
        cout << "Currently Frame: keyFrames[" << currentIndex << "].\n";
        frameToSgRbtNodes(*currentKeyFrame);
      }
      break;
    }


  case '+':
    if(g_msBetweenKeyFrames > 100){
      g_msBetweenKeyFrames-=100;
      cout << "accelerating..." << g_msBetweenKeyFrames << " ms between frames\n";
    }
    break;
  case '-': 
    if(g_msBetweenKeyFrames < 4000){
      g_msBetweenKeyFrames+=100;
      cout << "slowing..." << g_msBetweenKeyFrames << " ms between frames\n";
    }
    break;
  // - 'v': cycle through the 3 views
  case 'v':
    CHOOSEN_FRAME += (CHOOSEN_FRAME == 2) ? -2 : 1;
    cerr << "Currently active eye frame: " << MODE[CHOOSEN_FRAME] <<"\n";

    break;
  // - 'o': cycle through the 3 objects being manipulated

  case 'm':
    if(CHOOSEN_FRAME == 0){
      WORLD_SKY = !WORLD_SKY;
      cerr << "Switching to " << SKY_MODE[WORLD_SKY?0:1] << " frame" << "\n";
    }
    break;

  case 'p':
    PICKING = true;
    cout << "PICKING: " << "ON" << endl;
    return;

  // - 'r': reset the scene
  case 'r':
    reset();
    break;
  }
  glutPostRedisplay();
}


static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("Assignment 4");                       // title the window

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
}

static void initGLState() {
  glClearColor(128./255., 200./255., 255./255., 0.);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
  if (!g_Gl2Compatible)
    glEnable(GL_FRAMEBUFFER_SRGB);
}

static void initShaders() {
  g_shaderStates.resize(g_numShaders);
  for (int i = 0; i < g_numShaders; ++i) {
    if (g_Gl2Compatible)
      g_shaderStates[i].reset(new ShaderState(g_shaderFilesGl2[i][0], g_shaderFilesGl2[i][1]));
    else
      g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
  }
}

static void initGeometry() {
  initGround();
  initCubes();
  initSphere();
}

static void constructRobot(shared_ptr<SgTransformNode> base, const Cvec3& color) {

  const double ARM_LEN = 0.35,
               ARM_THICK = 0.25,
               L_ARM_LEN = 0.35,
               L_ARM_THICK = 0.125,
               HAND_LEN = 0.12,
               HEAD_LEN = 0.35,
               LEG_LEN = 0.5,
               LEG_THICK = 0.125,
               FOOT_LEN = 0.5,
               FOOT_THICK = 0.2,
               TORSO_LEN = 1.5,
               TORSO_THICK = 0.25,
               TORSO_WIDTH = 1;

  const int NUM_JOINTS = 12,
            NUM_SHAPES = 12;

  struct JointDesc {
    int parent;
    float x, y, z;
  };

  JointDesc jointDesc[NUM_JOINTS] = {
    {-1}, // torso
    {0, 0, TORSO_LEN/2, 0}, // head
    {0,  TORSO_WIDTH/2, TORSO_LEN/2, 0}, // upper right arm
    {0,  -TORSO_WIDTH/2, TORSO_LEN/2, 0}, // upper left arm
    {0, TORSO_WIDTH/4, -TORSO_LEN/2, 0}, // lower right leg
    {0, -TORSO_WIDTH/4, -TORSO_LEN/2, 0}, // lower left leg
    {2, ARM_LEN, 0, 0}, // upper right lower arm
    {3, -ARM_LEN, 0, 0}, // upper left lower arm
    {4, 0, -LEG_LEN*2, 0}, // lower right foot
    {5, 0, -LEG_LEN*2, 0}, // lower left foot
    {6, L_ARM_LEN*2, 0, 0}, // upper right hand
    {7, -L_ARM_LEN*2, 0, 0} // upper left hand
  };

  struct ShapeDesc {
    int parentJointId;
    float x, y, z, sx, sy, sz;
    shared_ptr<Geometry> geometry;
  };

  ShapeDesc shapeDesc[NUM_SHAPES] = {
    {0, 0, 0, 0, TORSO_WIDTH, TORSO_LEN, TORSO_THICK, g_cube}, // torso
    {1, 0, HEAD_LEN, 0, HEAD_LEN, HEAD_LEN, HEAD_LEN, g_sphere}, // head
    {2, ARM_LEN/2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube}, // upper right arm
    {3, -ARM_LEN/2, 0, 0, ARM_LEN, ARM_THICK, ARM_THICK, g_cube}, // upper left arm
    {4, 0, -LEG_LEN, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_sphere}, // lower right leg
    {5, 0, -LEG_LEN, 0, LEG_THICK, LEG_LEN, LEG_THICK, g_sphere}, // lower left leg
    {6, L_ARM_LEN, 0, 0, L_ARM_LEN, L_ARM_THICK, L_ARM_THICK, g_sphere}, // upper right lower arm
    {7, -L_ARM_LEN, 0, 0, L_ARM_LEN, L_ARM_THICK, L_ARM_THICK, g_sphere}, // upper left lower arm
    {8, 0, -FOOT_LEN/2, 0, FOOT_THICK, FOOT_LEN, FOOT_THICK, g_cube}, // lower right foot
    {9, 0, -FOOT_LEN/2, 0, FOOT_THICK, FOOT_LEN, FOOT_THICK, g_cube}, // lower left foot
    {10, HAND_LEN, 0, 0, HAND_LEN, HAND_LEN, HAND_LEN, g_sphere}, // upper right hand
    {11, -HAND_LEN, 0, 0, HAND_LEN, HAND_LEN, HAND_LEN, g_sphere} // upper left hand
  };

  shared_ptr<SgTransformNode> jointNodes[NUM_JOINTS];

  for (int i = 0; i < NUM_JOINTS; ++i) {
    if (jointDesc[i].parent == -1)
      jointNodes[i] = base;
    else {
      jointNodes[i].reset(new SgRbtNode(RigTForm(Cvec3(jointDesc[i].x, jointDesc[i].y, jointDesc[i].z))));
      jointNodes[jointDesc[i].parent]->addChild(jointNodes[i]);
    }
  }
  for (int i = 0; i < NUM_SHAPES; ++i) {
    shared_ptr<MyShapeNode> shape(
      new MyShapeNode(shapeDesc[i].geometry,
                      color,
                      Cvec3(shapeDesc[i].x, shapeDesc[i].y, shapeDesc[i].z),
                      Cvec3(0, 0, 0),
                      Cvec3(shapeDesc[i].sx, shapeDesc[i].sy, shapeDesc[i].sz)));
    jointNodes[shapeDesc[i].parentJointId]->addChild(shape);
  }
}

static void initScene() {
  g_world.reset(new SgRootNode());

  g_skyNode.reset(new SgRbtNode(RigTForm(Cvec3(0.0, 0.25, 4.0))));

  g_groundNode.reset(new SgRbtNode());
  g_groundNode->addChild(shared_ptr<MyShapeNode>(
                           new MyShapeNode(g_ground, Cvec3(0.1, 0.95, 0.1))));

  g_robot1Node.reset(new SgRbtNode(RigTForm(Cvec3(-2, 1, 0))));
  g_robot2Node.reset(new SgRbtNode(RigTForm(Cvec3(2, 1, 0))));

  constructRobot(g_robot1Node, Cvec3(1, 0, 0)); // a Red robot
  constructRobot(g_robot2Node, Cvec3(0, 0, 1)); // a Blue robot

  g_world->addChild(g_skyNode);
  g_world->addChild(g_groundNode);
  g_world->addChild(g_robot1Node);
  g_world->addChild(g_robot2Node);

  // dump all the nodes to nodeData
  dumpSgRbtNodes(g_world->shared_from_this(), nodeData);
  // current picked node is NULL
  g_currentPickedRbtNode = shared_ptr<SgRbtNode>(); 

}

int main(int argc, char * argv[]) {
  try {
    initGlutState(argc,argv);

    glewInit(); // load the OpenGL extensions

    cout << (g_Gl2Compatible ? "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.3") << endl;
    if ((!g_Gl2Compatible) && !GLEW_VERSION_3_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
    else if (g_Gl2Compatible && !GLEW_VERSION_2_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.0");

    initGLState();
    initShaders();
    initGeometry();
    initScene();

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}

// void printRigTForm(RigTForm p){
//   cout << "Cvec3: ";
//   for(int i = 0;i < 3;i++){
//     cout << p.getTranslation()[i] << " ";
//   }
//   cout << "|| Quat4: ";
//   for(int i = 0;i < 4;i++){
//     cout << p.getRotation()[i] << " ";
//   }
//   cout << endl;
// }

