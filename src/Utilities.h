/***********************************************************
*   Simple GLSL demo code for COS426 Computer Graphics     *
*                           ****                           *
*   Most is borrowed from                                  *
*   http://www.lighthouse3d.com/opengl/glsl/               *
*                           ****                           *
*   Dependencies: glew, glut, gl                           *
*                           ****                           *
*   Pulled together by Aleksey Boyko(aboyko@princeton.edu) *
************************************************************/

#ifndef __GLSL_TUTORIAL_UTILITIES__
#define __GLSL_TUTORIAL_UTILITIES__

#include <iostream>
#include <GL/glew.h>
#if defined(_WIN32) || defined(__CYGWIN__)
# include <windows.h>
# include <GL/glut.h>
#elif defined(__APPLE__)
# include <GLUT/glut.h>
#else 
# include <GL/glut.h>
#endif

using namespace std;

float initGlew(bool verbose = false);

void printShaderInfoLog(GLuint obj);
void printProgramInfoLog(GLuint obj);

char *textFileRead(const char *fn);
int textFileWrite(const char *fn, const char *s);

#endif

