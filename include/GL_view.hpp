#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include <iostream>

#include <boost/timer/timer.hpp>

/* OpenGL and friends */
#ifdef USEGLEW
#include <GL/glew.h>
#endif
#define GL_GLEXT_PROTOTYPES

#include <GL/glut.h>
#include "Matrix4D.hpp"
/*  Macro for sin & cos in degrees */
#define Cos(th) cos(M_PI/180*(th))
#define Sin(th) sin(M_PI/180*(th))

#define LEN 8192

static double dim = 8;
static char *windowName = "Rotation, Scaling, and Translations";
static int windowWidth = 1280;
static int windowHeight = 900;

static int toggleAxes = 1;
static int toggleValues = 1;
static int toggleMode = 0;
static std::vector<float> error_{};
static int th = 0;
static int ph = 0;
static int fov = 55;
static int asp = 1;
static float t{};

static constexpr int sampleMin = 0;
static constexpr int sampleMax = 5;

static constexpr int NI = 6;
static constexpr int NJ = 6;
static constexpr int NK = 6;
static constexpr int NL = 6;

using namespace _spatial;
typedef boost::multi_array<int, 4> array_type;
static const std::array<array_type::index, 4> shape = {{NI, NJ, NK, NL}};
static constexpr std::array<array_type::index, 4> shape1{4, 5, 2, 3};
static constexpr std::array<array_type::index, 4> shape2{4, 5, 3, 4};

inline float deg2rad (float deg)
{
    return deg*M_PI/180;
}

inline void printv(va_list args, const char* format)
{
  char buf[LEN];
  char* ch=buf;
  vsnprintf(buf,LEN,format,args);
  while (*ch)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*ch++);
}

inline void print(const char* format, ...)
{
  va_list args;
  va_start(args,format);
  printv(args,format);
  va_end(args);
}

inline void printAt(int x, int y, const char* format, ...)
{
  va_list args;
  glWindowPos2i(x,y);
  va_start(args,format);
  printv(args,format);
  va_end(args);
}


inline void errCheck(char* where)
{
   int err = glGetError();
   if (err) fprintf(stderr,"ERROR: %s [%s]\n",gluErrorString(err),where);
}


inline void glWrite(float x, float y, int *font, const char* text, size_t kls) {
    glRasterPos2f(x, y);
    for (size_t i = 0; i < kls; i++)
    glutBitmapCharacter(font, text[i]);
}

inline void project()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (toggleMode) {
    /* perspective */
    gluPerspective(fov, asp, dim/4, 4*dim);
  }
  else {
    /* orthogonal projection*/
    glOrtho(-dim*asp, +dim*asp, -dim, +dim, -dim, +dim);
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

inline void setEye()
{
  if (toggleMode) {
    double Ex = -2*dim*Sin(th)*Cos(ph);
    double Ey = +2*dim        *Sin(ph);
    double Ez = +2*dim*Cos(th)*Cos(ph);

    gluLookAt(Ex,Ey,Ez , 0,0,0 , 0,Cos(ph),0);
  }  else {
    glRotatef(ph,1,0,0);
    glRotatef(th,0,1,0);
  }
}


inline void drawValues()
{
  if (toggleValues) {
    glColor3f(0.8,0.8,0.8);
    printAt(5,5,"View Angle (th, ph) =(%d, %d)", th,ph);
    printAt(5,25,"Projection mode =(%s)", toggleMode?"Perspective":"Orthogonal");
        for(size_t i = 0; i < error_.size(); ++i)
    printAt(5, 880 - 20*i,"measurement error |f(x, y) - p(x, y)| =(%f)", error_[i]);
    error_.clear();
  }
}


inline void initlights(void)
{
   GLfloat ambient[] = {0.6, 0.6, 0.6, 1.0};
   GLfloat position[] = {-16.0, 16.0, 10.0, .3};
   GLfloat mat_diffuse[] = {.8, .8, .8, 1.0};
   GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
   GLfloat mat_shininess[] = {80.0};

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);

   glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
   glLightfv(GL_LIGHT0, GL_POSITION, position);

   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
}


inline void drawSurface() {
    try {
        boost::timer::cpu_timer tmr;
        Matrix4D<int> A(shape), B(shape), C(shape);

        A.Random(1, 5);
        B.Random(2, 6);

        C +=  A + (B+5);

        glPointSize(10);
        auto ptr_trns = std::make_unique<decltype (C.transversal_section('i'))>(C.transversal_section('i'));

        glBegin (GL_POINTS);
            for(const auto& x1 : ptr_trns->second) {
                glColor3f(1, 1, 0); glVertex3f (x1[0] + x1[3], x1[1] + x1[3], x1[2] + x1[3]);
            }
        glEnd ();
        glBegin(GL_LINES);
              for(size_t i = 0; i < A.size(0); i += 1 ) {
                  for(size_t j = 0; j < A.size(1); j += 1 ) {
                      for (size_t k = 0; k < A.size(2); k ++) {
                          for (size_t l = 0; l < A.size(3); l ++) {
                              glColor3f(1, 0, 0); glVertex3f(l, j + l, k + l); glVertex3f(A.size(0) - 1 + l, j + l, k + l);
                              glColor3f(1, 0, 0); glVertex3f(i + l, l, k + l); glVertex3f(i + l, A.size(1) -1 + l, k + l);
                              glColor3f(1, 0, 0); glVertex3f(i + l, j + l, l); glVertex3f(i + l, j + l, A.size(2) - 1 + l);
                          }
                      }

                  }
              }
          glEnd();

          for(const auto& x1 : ptr_trns->second) {
            glColor3f(1.0, 1.0, 1.0);
            glRasterPos3d(x1[0] + x1[3], x1[1] + x1[3], x1[2] + x1[3]);
            char buffer[80];
            sprintf(buffer, "(%ld,%ld,%ld,%d)",x1[0], x1[1], x1[2], x1[3]/*, C(i, j, k)*/);
            print(buffer);
          }
          glFlush();
          std::cout << "UPDATE: " << tmr.format() << std::endl;

    }  catch (std::exception& e) {
        std::cerr << "ERROR_poison***********" << e.what() << std::endl;
    }

}

inline void display()
{
    glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    setEye();
    drawValues();

    drawSurface();
    glPopMatrix();
    glFlush();
    glutSwapBuffers();
}


inline void reshape(int width,int height)
{
    asp = (height>0) ? (double)width/height : 1;
    glViewport(0,0, width,height);
    project();

}


inline void windowKey(unsigned char key, int x, int y)
{
    /*  Exit on ESC */
    if (key == 27) exit(0);
    else if (key == 'a' || key == 'A') toggleAxes = 1-toggleAxes;
    else if (key == 'v' || key == 'V') toggleValues = 1-toggleValues;
    else if (key == 'm' || key == 'M') toggleMode = 1-toggleMode;
    /*  Change field of view angle */
    else if (key == '-' && key>1) fov--;
    else if (key == '+' && key<179) fov++;
    /*  Change dimensions */
    else if (key == 'D') dim += 0.1;
    else if (key == 'd' && dim>1) dim -= 0.1;

    project();
    glutPostRedisplay();
}


inline void windowSpecial(int key,int x,int y)
{
  /*  Right arrow key - increase azimuth by 5 degrees */
  if (key == GLUT_KEY_RIGHT) th += 5;
  /*  Left arrow key - decrease azimuth by 5 degrees */
  else if (key == GLUT_KEY_LEFT) th -= 5;
  /*  Up arrow key - increase elevation by 5 degrees */
  else if (key == GLUT_KEY_UP) ph += 5;
  /*  Down arrow key - decrease elevation by 5 degrees */
  else if (key == GLUT_KEY_DOWN) ph -= 5;

  /*  Keep angles to +/-360 degrees */
  th %= 360;
  ph %= 360;

  project();
  glutPostRedisplay();
}

inline void MouseWheel(int wheel, int direction, int x, int y)
{
    wheel=0;
    if (direction==-1)
    {
        dim -= 0.5;

    }
    else if (direction==+1)
    {
        dim += 0.5;
    }
project();
 glutPostRedisplay();

}
