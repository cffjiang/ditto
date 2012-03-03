//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/visualization/viewer.h
// Copyright 2012, Chenfanfu Jiang
//
// OpenGL viewer.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_VISUALIZATION_VIEWER_H
#define DITTO_PUBLIC_LIBRARY_VISUALIZATION_VIEWER_H

#include <GL/glut.h>
#include <cstdlib>
#include <cmath>
#include "arcball.h"

namespace ditto{ namespace visualization {

GLint uAmbient, uDiffuse, uSpecular, uLightPos, uShininess;

static float aspect_ratio = 1.0f;
static int width, height;

// scene parameters
vec eye( 0.0f, 0.0f, -20.0f );
vec centre( 0.0f, 0.0f, 0.0f );
vec object_translate( 0.0f, 0.0f, 0.0f );
const vec up( 0.0f, 1.0f, 0.0f );
const float SPHERE_RADIUS = 5.0f;
const int SPHERE_LAT_SLICES = 12;
const int SPHERE_LONG_SLICES = 24;

const float PI = 3.141592654f;

static void reset_view(int w, int h);
inline float randf();
static void shutdown_scene();
inline vec rotate_x( vec v, float sin_ang, float cos_ang );
inline vec rotate_y( vec v, float sin_ang, float cos_ang );
static void draw_sphere();
static void click_scene(int x, int y);
static void drag_scene(int x, int y);
static void draw_scene();
static void resize(int w, int h);
static void display();
static void key(unsigned char key, int x, int y);
static void idle();
static void mouse_button(int button, int state, int x, int y);
static void mouse_motion(int x, int y);
int main_loop(int argc, char ** argv);

static void reset_view(int w, int h)
{
	width = w;
	height = h;
    aspect_ratio = (float) width / (float) height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( 50.0f, aspect_ratio, 1.0f, 50.0f );
    gluLookAt(
        eye.x, eye.y, eye.z,
        centre.x, centre.y, centre.z,
        up.x, up.y, up.z );
	// set up the arcball using the current projection matrix
	arcball_setzoom( SPHERE_RADIUS, eye, up );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;
}

void set_colour(float r, float g, float b)
{
    float ambient = 0.2f;
    float diffuse = 0.7f;
    float specular = 0.7f;
    GLfloat mat[4];
    /**** set ambient lighting parameters ****/
    mat[0] = ambient*r;
    mat[1] = ambient*g;
    mat[2] = ambient*b;
    mat[3] = 1.0;
    glMaterialfv (GL_FRONT, GL_AMBIENT, mat);

    /**** set diffuse lighting parameters ******/
    mat[0] = diffuse*r;
    mat[1] = diffuse*g;
    mat[2] = diffuse*b;
    mat[3] = 1.0;
    glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);

    /**** set specular lighting parameters *****/
    mat[0] = specular*r;
    mat[1] = specular*g;
    mat[2] = specular*b;
    mat[3] = 1.0;
    glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
    glMaterialf (GL_FRONT, GL_SHININESS, 1.0);
}

void my_init()
{
    // GLfloat ambient[] = { 0.2, 0.2, 0.2, 1.0 };
    GLfloat ambient[] = { 0.8, 0.8, 0.8, 1.0 };
    GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
    // GLfloat position[] = { -3.0, 3.0, 3.0, 0.0 };
    GLfloat position[] = { 10.0, 0.0, -30.0, 1.0 };
    GLfloat diffuse2[] = { 0.3, 0.3, 0.3, 1.0 };
    GLfloat specular2[] = { 0.3, 0.3, 0.3, 1.0 };
    GLfloat position2[] = { 0.0, 100.0, 00.0, 1.0 };

    GLfloat lmodel_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat local_view[] = { 0.0 };

    /**** set lighting parameters ****/
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    glLightf(GL_LIGHT0, GL_SHININESS, 100) ;
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE) ;

    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse2);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular2);
    glLightfv(GL_LIGHT1, GL_POSITION, position2);
    glLightf(GL_LIGHT1, GL_SHININESS, 500) ;

    /*    glFrontFace (GL_CW); */
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glPixelStorei(GL_PACK_ALIGNMENT,1) ;
    glPixelStorei(GL_UNPACK_ALIGNMENT,1) ;
    glShadeModel(GL_SMOOTH) ;
}

inline float randf()
{
	return ((1.0f / 127.f) * (((float)(rand() & 255)) - 128.0f)) ;
}

static void shutdown_scene()
{
	// nothing to be done here
}

inline vec rotate_x( vec v, float sin_ang, float cos_ang )
{
	return vec(
	    v.x,
	    (v.y * cos_ang) + (v.z * sin_ang),
	    (v.z * cos_ang) - (v.y * sin_ang)
	    );
}

inline vec rotate_y( vec v, float sin_ang, float cos_ang )
{
	return vec(
	    (v.x * cos_ang) + (v.z * sin_ang),
	    v.y,
	    (v.z * cos_ang) - (v.x * sin_ang)
	    );
}

static void draw_sphere()
{
    glutSolidSphere(1,30,20);
}

static void click_scene(int x, int y)
{
	int invert_y = (height - y) - 1; // OpenGL viewport coordinates are Cartesian
	arcball_start(x,invert_y);
}

static void drag_scene(int x, int y)
{
	int invert_y = (height - y) - 1;
	arcball_move(x,invert_y);
}

static void draw_scene()
{
    glTranslatef( object_translate.x, object_translate.y, object_translate.z );
    arcball_rotate();
    
    // draw things below
    
    set_colour( 0.0f, 0.0f, 1.0f );
	draw_sphere();
}

static void resize(int w, int h)
{
	reset_view(w,h);
}

static void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    draw_scene();

    glutSwapBuffers();
}

static void key(unsigned char key, int x, int y)
{
    switch (key) 
    {
        case 27 : 
        case 'q':
			shutdown_scene();
            exit(0);
            break;
        case 'w':
            object_translate.y += .2;
            break;
        case 's':
            object_translate.y -= .2;
            break;
        case 'a':
            object_translate.x += .2;
            break;
        case 'd':
            object_translate.x -= .2;
            break;
        case '-':
            object_translate.z += .2;
            break;
        case '+':
            object_translate.z -= .2;
            break;
		default :
			break;
    }
    glutPostRedisplay();
}

static void idle()
{
    glutPostRedisplay();
}

static void mouse_button(int button, int state, int x, int y)
{
	if ( state == GLUT_DOWN ) click_scene(x,y);
}

static void mouse_motion(int x, int y)
{
	// glutMotionFunc only called when a mouse button is held down
	drag_scene(x,y);
}

int main_loop(int argc, char ** argv)
{
	arcball_reset();

    glutInit(&argc, argv);
    glutInitWindowSize(800,600);
    glutInitWindowPosition(10,10);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Ditto OpenGL Viewer");

    my_init();
    
    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutIdleFunc(idle);
    glutMouseFunc(mouse_button);
    glutMotionFunc(mouse_motion);

    glClearColor(0.0f,0.0f,0.0f,1.0f);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glutMainLoop();

    return 0;
}

} }

#endif
