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

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
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
static const bool g_Gl2Compatible = true;


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
// ========================================
// TODO: you can add global variables here
// ========================================

// hw2 globals -------------------------------------------------
// number of cubes to draw
static const int CUBE_NUMBER = 2;
// default, choose sky camera
// sky camera: 0; cube one: 1; cube two: 2;
static int CHOOSEN_OBJECT = 0;
// default, sky camera frame
// sky camera frame: 0; cube one: 1; cube two: 2;
static int CHOOSEN_FRAME = 0;
// default, world sky camera
static bool WORLD_SKY = true;
// just names
static string MODE[] = {"sky", "cube 1", "cube 2"};
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


struct ShaderState {
  GlProgram program;

  // Handles to uniform variables
  GLint h_uLight, h_uLight2;
  GLint h_uProjMatrix;
  GLint h_uModelViewMatrix;
  GLint h_uNormalMatrix;
  GLint h_uColor;

  // Handles to vertex attributes
  GLint h_aPosition;
  GLint h_aNormal;

  ShaderState(const char* vsfn, const char* fsfn) {
    readAndCompileShader(program, vsfn, fsfn);

    const GLuint h = program; // short hand

    // Retrieve handles to uniform variables
    h_uLight = safe_glGetUniformLocation(h, "uLight");
    h_uLight2 = safe_glGetUniformLocation(h, "uLight2");
    h_uProjMatrix = safe_glGetUniformLocation(h, "uProjMatrix");
    h_uModelViewMatrix = safe_glGetUniformLocation(h, "uModelViewMatrix");
    h_uNormalMatrix = safe_glGetUniformLocation(h, "uNormalMatrix");
    h_uColor = safe_glGetUniformLocation(h, "uColor");

    // Retrieve handles to vertex attributes
    h_aPosition = safe_glGetAttribLocation(h, "aPosition");
    h_aNormal = safe_glGetAttribLocation(h, "aNormal");

    if (!g_Gl2Compatible)
      glBindFragDataLocation(h, 0, "fragColor");
    checkGlErrors();
  }

};

static const int g_numShaders = 2;
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"}
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
  {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"}
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


// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube;
// hw3 a shared sphere
static shared_ptr<Geometry> g_sphere;

// --------- Scene

static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
//static Matrix4 g_skyRbt = Matrix4::makeTranslation(Cvec3(0.0, 0.25, 4.0)); // hw2 Matrix4 g_skyRbt
static RigTForm g_skyRbt(Cvec3(0.0, 0.25, 4.0)); // hw3 RigTForm g_skyRbt
// ============================================
// TODO: add a second cube's 
// 1. transformation
// 2. color
// ============================================

// hw2 globals --------------------------------
// Added a second cube 
// Two cube's coordinates are 1st (-1,0,0) and 2nd (1,0,0)
// Two cube's colors are 1st pure red(on the left) and 2nd pure blue(on the right)
// --------------------------------------------

// static Matrix4 g_objectRbt[CUBE_NUMBER] = {Matrix4::makeTranslation(Cvec3(-1,0,0)), 
//                                            Matrix4::makeTranslation(Cvec3(1,0,0))}; 

static RigTForm g_objectRbt[CUBE_NUMBER] = {RigTForm(Cvec3(-1,0,0)), 
                                            RigTForm(Cvec3(1,0,0))};

static Cvec3f g_objectColors[CUBE_NUMBER] = {Cvec3f(1, 0, 0), 
                                             Cvec3f(0, 0, 1)};

// hw3 global
// sphere center
static Cvec3 sphere_center_eyeCord; 


///////////////// END OF G L O B A L S //////////////////////////////////////////////////




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

// takes MVM and its normal matrix to the shaders
static void sendModelViewNormalMatrix(const ShaderState& curSS, const Matrix4& MVM, const Matrix4& NMVM) {
  GLfloat glmatrix[16];
  MVM.writeToColumnMajorMatrix(glmatrix); // send MVM
  safe_glUniformMatrix4fv(curSS.h_uModelViewMatrix, glmatrix);

  NMVM.writeToColumnMajorMatrix(glmatrix); // send NMVM
  safe_glUniformMatrix4fv(curSS.h_uNormalMatrix, glmatrix);
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

static void drawStuff() {
  // short hand for current shader state
  const ShaderState& curSS = *g_shaderStates[g_activeShader];

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  // Choose an eye camera
  //Matrix4 cam;
  RigTForm cam;
  switch(CHOOSEN_FRAME){
    case 0:
      cam = g_skyRbt;
      break;
    case 1:
      cam = g_objectRbt[0];
      break;
    case 2:
      cam = g_objectRbt[1];
      break;
  }

  //const Matrix4 eyeRbt = cam;
  //const Matrix4 invEyeRbt = inv(eyeRbt);
  const RigTForm eyeRbt = cam;
  const RigTForm invEyeRbt = inv(eyeRbt);

  const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
  safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
  safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

  // draw ground
  // ===========
  //
  //const Matrix4 groundRbt = Matrix4();  // identity
  const RigTForm groundRbt = RigTForm();  // identity
  //Matrix4 MVM = invEyeRbt * groundRbt;
  Matrix4 MVM = rigTFormToMatrix(invEyeRbt * groundRbt);
  Matrix4 NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, 0.1, 0.95, 0.1); // set color
  g_ground->draw(curSS);

  // draw cubes (now two)
  // ----------
  for (int i = 0;i < CUBE_NUMBER;i++){
    //MVM = invEyeRbt * g_objectRbt[i];
    MVM = rigTFormToMatrix(invEyeRbt * g_objectRbt[i]);
    NMVM = normalMatrix(MVM);
    sendModelViewNormalMatrix(curSS, MVM, NMVM);
    safe_glUniform3f(curSS.h_uColor, g_objectColors[i][0], g_objectColors[i][1], g_objectColors[i][2]);
    g_cube->draw(curSS);
  }

  // draw a sphere?
  if(!DRAW_SPHERE){return;}

  // when to draw a sphere
  // 1. when in world sky system manipulating sky camera
  // 2. when manipulating a cube but not in its own frame
  const bool camera_world_sky = !CHOOSEN_FRAME && !CHOOSEN_OBJECT && WORLD_SKY;
  const bool cube_not_itself = CHOOSEN_OBJECT && (CHOOSEN_FRAME != CHOOSEN_OBJECT);

  if(camera_world_sky || cube_not_itself){

    RigTForm sphere_center;// now identity
    switch(CHOOSEN_OBJECT){
      case 0:
        break;// if in world sky manipulating sky camera, center should be in world center, a.k.a identity()
      case 1:
        sphere_center = g_objectRbt[0];
        break;
      case 2:
        sphere_center = g_objectRbt[1];
        break;
    }

    // update the radius
    RigTForm eyeCord = invEyeRbt * sphere_center; 
    // refresh sphere center cordinates
    sphere_center_eyeCord = eyeCord.getTranslation();
    // ---
    // snap-back-to-constant-size effect
    if (!(g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton))){ 
      g_arcballScale = getScreenToEyeScale(sphere_center_eyeCord[2], g_frustFovY, g_windowHeight);
    }
    double arcball_radius = g_arcballScale * g_arcballScreenRadius;

    // draw sphere
    // ----------
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    MVM = rigTFormToMatrix(eyeCord) * Matrix4::makeScale(Cvec3(arcball_radius, arcball_radius, arcball_radius));
    NMVM = normalMatrix(MVM);
    sendModelViewNormalMatrix(curSS, MVM, NMVM);
    safe_glUniform3f(curSS.h_uColor, 0, 0, 0);
    g_sphere->draw(curSS);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  }
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff();

  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

  checkGlErrors();
}

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

static void motion(const int x, const int y) {
  const double dx = x - g_mouseClickX;
  const double dy = g_windowHeight - y - 1 - g_mouseClickY;

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

  // if there is no move, just return
  if(v0[0] == v1[0] && v0[1] == v1[1] && v0[2] == v1[2]){return;}

  // cross product
  // Cvec3 k = normalize(cross(v0, v1));
  // // because it is arcball, we spin twice the angle!
  // double angle_cos = dot(v0, v1) / (norm(v0)*norm(v1));
  // double angle_sin = sqrt(1 - angle_cos*angle_cos);

  //Quat theQuat(angle_cos, k*angle_sin);

  Quat theQuat(normalize(Quat(0, v1)) * normalize(Quat(0, -v0)));


  // Matrix4 Q, m, A;
  // Matrix4 objectFrame;
  // Matrix4 eyeFrame;

  RigTForm Q, A;
  RigTForm objectFrame;
  RigTForm eyeFrame;

  // choose transformation matrix and rotation matrix based on object we are manipulating and frame
  if((CHOOSEN_OBJECT == 1 && CHOOSEN_FRAME == 0) || (CHOOSEN_OBJECT == 2 && CHOOSEN_FRAME == 0)){
    objectFrame = g_objectRbt[CHOOSEN_OBJECT - 1];
    eyeFrame = g_skyRbt;
  }
  else if((CHOOSEN_OBJECT >= 1 && CHOOSEN_OBJECT <= 2) && (CHOOSEN_FRAME >= 1 || CHOOSEN_FRAME <= 2)){
    objectFrame = g_objectRbt[CHOOSEN_OBJECT - 1];
    eyeFrame = g_objectRbt[CHOOSEN_FRAME - 1];
  }
  else if(CHOOSEN_OBJECT == 0 && CHOOSEN_FRAME == 0){
    if(!WORLD_SKY)
      objectFrame = g_skyRbt;
    eyeFrame = g_skyRbt;
  }
  // calculate the auxiliary matrix
  A = makeMixedFrame(objectFrame, eyeFrame);
  // calculate the translation or rotation matrixs
  if (g_mouseLClickButton && !g_mouseRClickButton) { // left button down?
    //Q = Matrix4::makeXRotation(-dy) * Matrix4::makeYRotation(dx);
    //Q = RigTForm(Quat::makeXRotation(-dy) * Quat::makeYRotation(dx));
    Q = RigTForm(theQuat);
    if(CHOOSEN_OBJECT != 0 && CHOOSEN_FRAME == CHOOSEN_OBJECT)
      Q = inv(Q);
  }
  else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
    if(DRAW_SPHERE)
      Q = RigTForm(Cvec3(dx, dy, 0) * g_arcballScale);//hw3
    else
      Q = RigTForm(Cvec3(dx, dy, 0) * 0.01); 
    if(!WORLD_SKY && !CHOOSEN_FRAME && !CHOOSEN_OBJECT)
      Q = inv(Q);
  }
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
    if(DRAW_SPHERE)
      Q = RigTForm(Cvec3(0, 0, -dy) * g_arcballScale);//hw3
    else
      Q = RigTForm(Cvec3(0, 0, -dy) * 0.01);
    if(!WORLD_SKY && !CHOOSEN_FRAME && !CHOOSEN_OBJECT)
      Q = inv(Q);
  }

  if (g_mouseClickDown) {
    // DO THIS: O = A * Q * A^-1 * o
    switch(CHOOSEN_OBJECT){
    case 0:
      g_skyRbt = doQtoOwrtA(inv(Q), g_skyRbt, A);
      break;
    case 1:
      g_objectRbt[0] = doQtoOwrtA(Q, g_objectRbt[0], A);
      break;
    case 2:
      g_objectRbt[1] = doQtoOwrtA(Q, g_objectRbt[1], A);
      break;
    }
    
    glutPostRedisplay(); // we always redraw if we changed the scene
  }

  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}

static void reset()
{
  // =========================================================
  // TODO:
  // - reset g_skyRbt and g_objectRbt to their default values
  // - reset the views and manipulation mode to default
  // - reset sky camera mode to use the "world-sky" frame
  // =========================================================

  //reset g_skyRbt and g_objectRbt to their default values
  g_skyRbt = RigTForm(Cvec3(0.0, 0.25, 4.0));
  g_objectRbt[0] = RigTForm(Cvec3(-1,0,0));
  g_objectRbt[1] = RigTForm(Cvec3(1,0,0));
  //reset the views and manipulation mode to default
  CHOOSEN_OBJECT = 0;
  CHOOSEN_FRAME = 0;
  //reset sky camera mode to use the "world-sky" frame
  WORLD_SKY = true;
  // draw sphere
  DRAW_SPHERE = true;

	cout << "reset objects and modes to defaults" << endl;
  cout << "reset to draw sphere mode" << endl;
}

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
  glutPostRedisplay();
}


static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC

  // -'a' control whether draw the sphere
  case 'a':
    DRAW_SPHERE = !DRAW_SPHERE;

    cout << "Now in ";
    if(DRAW_SPHERE)cout << "SPHERE";
    else cout << "NO-SPHERE";  
    cout << " mode.\n";
    
    break;
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "o\t\tCycle object to edit\n"
    << "v\t\tCycle view\n"
    << "drag left mouse to rotate\n" << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;
  case 'f':
    g_activeShader ^= 1;
    break;
  // ============================================================
  // TODO: add the following functionality for 
  //       keybaord inputs
  // - 'v': cycle through the 3 views
  // - 'o': cycle through the 3 objects being manipulated
  // - 'm': switch between "world-sky" frame and "sky-sky" frame
  // - 'r': reset the scene
  // ============================================================

  // - 'v': cycle through the 3 views
  case 'v':
    CHOOSEN_FRAME += (CHOOSEN_FRAME == 2) ? -2 : 1;
    // When changing a frame, change the object to frame itself
    // it will be easier to operate
    CHOOSEN_OBJECT = CHOOSEN_FRAME;
    cerr << "Currently active eye frame: " << MODE[CHOOSEN_FRAME] <<"\n";
    break;
  // - 'o': cycle through the 3 objects being manipulated
  case 'o':
    if(CHOOSEN_FRAME == 0)
      CHOOSEN_OBJECT += (CHOOSEN_OBJECT == 2) ? -2 : 1;
    // prevent changing to camera object while in the cube view
    else
      CHOOSEN_OBJECT += (CHOOSEN_OBJECT == 2) ? -1 : 1;
    break;
  // - 'm': switch between "world-sky" frame and "sky-sky" frame
  case 'm':
    if(CHOOSEN_FRAME == 0 && CHOOSEN_OBJECT == 0){
      WORLD_SKY = !WORLD_SKY;
      cerr << "Switching to " << SKY_MODE[WORLD_SKY?0:1] << " frame" << "\n";
    }
    break;
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
  glutCreateWindow("Assignment 2");                       // title the window

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

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
