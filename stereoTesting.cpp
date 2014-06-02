/*
 *
 *  Filename:  stereoTesting.cpp
 *
 *  Description:
 *    This code is to demonstrate how we have implemented the stereoscopic
 *    viewing for the visweek paper. 
 *
 *    The code is written on a Linux platform and the computer we used to
 *    run the study.
 *
 * Creation: May 25 2014
 *
 * Author: Jian Chen
 * 
 *  Usage:
 *    ./stereoTest 
 *        will display mono 
 *    ./stereoTest -stereo -f ./data/test.data
 *
 *  Requirement: 
 *    - test.data: a brain tractography dataset.
 *    - your graphics card must support stereoscopic display. The graphics
 *      card we used in the experiment was a nVidia Quadro 4000 (by PNY). 
 *    - you will also need stereoscopic glasses, e.g., the 3D vision kit
 *      also from nVidia.
 *
 * Acknowledgement:
 *    This method is modified from the stereo rendering code written by 
 *    Paul Bourke. We have adapted it to our experimental conditions
 *    in showing the brain dataset.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>
#include "GLoader.h"

#define PI		M_PI
#define DTOR            (M_PI/180.0)
#define RTOD            (180.0 / M_PI) //57.2957795
#define CROSSPROD(p1,p2,p3) \
   p3.x = p1.y*p2.z - p1.z*p2.y; \
   p3.y = p1.z*p2.x - p1.x*p2.z; \
   p3.z = p1.x*p2.y - p1.y*p2.x

// bounding volume of the brain model
GLfloat Xminmax[2] = {50.774986, 182.695862};
GLfloat Yminmax[2] = {38.447567, 203.434158};
GLfloat Zminmax[2] = {5.271503, 135.014359};
CStreamtubes tubes;

int datafile = GL_FALSE;

GLint rot_x = 0, rot_y = 0, rot_z = 0;
GLfloat scale= 3.;

typedef struct {
    double x, y, z;
} XYZ;
typedef struct {
    double r, g, b;
} COLOUR;
typedef struct {
    unsigned char r, g, b, a;
} PIXELA;
typedef struct {
    XYZ vp;			/* View position           */
    XYZ vd;			/* View direction vector   */
    XYZ vu;			/* View up direction       */
    XYZ pr;			/* Point to rotate about   */
    double focallength;		/* Focal Length along vd   */
    double aperture;		/* Camera aperture         */
    double eyesep;		/* Eye separation          */
    int screenwidth, screenheight;
} CAMERA;
double dtheta = 1;
CAMERA camera;
XYZ origin = { 0.0, 0.0, 0.0 };

/*------------------------------------------------------------------------
 * Normalize:
 *   Input:
 *   Output:
 *   Description:
 *      handy math routine to calculate vector normal.
 */
static void Normalise(XYZ * p)
{
    double length;

    length = sqrt(p->x * p->x + p->y * p->y + p->z * p->z);
    if (length != 0) {
		p->x /= length;
		p->y /= length;
		p->z /= length;
    } else {
		p->x = 0;
		p->y = 0;
		p->z = 0;
    }
}

/*------------------------------------------------------------------------
 *  initCamera
 *    Input:
 *    Output:
 *    Description:
 *       camera setup
 *       FOV = 60;
 *       Focal_length = viewing distance = 460.1mm on the small 24" display
 *  
 */
static void initCamera()
{
    camera.aperture = 60;
    camera.focallength = 460.1;
    camera.eyesep = camera.focallength / 30;  // for comfortable viewing
                                              // the denominator can be
                                              // any number between 20 and 30.
    camera.pr = origin;

    camera.vp.x = 0;
    camera.vp.y = 0;
    camera.vp.z = (Zminmax[0]+Zminmax[1]) * 6.0;
    camera.vd.x = -camera.vp.x;
    camera.vd.y = -camera.vp.y;
    camera.vd.z = -camera.vp.z;

    camera.vu.x = 0;
    camera.vu.y = 1;
    camera.vu.z = 0;


    camera.screenwidth = 1920;
    camera.screenheight = 1080;
}

/*------------------------------------------------------------------------
 * KeyFunc ()
 *  Input:
 *  Output:
 *  Description:
 *  	Some handy keyinput.
 */
static void keyFunc(unsigned char key, int x, int y)
{
    switch (key) {
		case 27:			/* Quit */
		case 'Q':
		case 'q':
			exit(0);
		case 'x':
			rot_x = (rot_x + 1) % 360;
			break;
		case 'y':
			rot_y = (rot_y + 1) % 360;
			break;
		case 'z':
			rot_z = (rot_z + 1) % 360;
			break;
		case 32: // SPACE to get home
			rot_x = rot_y = rot_z = 0;
			break;
		case 'c':
			glDumpImagePPM(0,0,1920,1080,true);
			break;
		case '=':
			scale += 0.1;
			break;
		case '-':
			scale -= 0.1;
			break;
		default:
			return;
	}
	glutPostRedisplay();
}

/*------------------------------------------------------------------------
 * do_draw()
 *   Input: none
 *   Output: 
 *      geometry drawings. See draw() for details
 */
static void do_draw(void)
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glClearColor(0.,0.,0.,0);

   glScalef(scale, scale, scale);

   // three axes: rotate with the model 
   glPushMatrix();
     glRotatef(rot_x, 1,0,0);
     glRotatef(rot_y, 0,1,0);
     glRotatef(rot_z, 0,0,1);

     glBegin(GL_LINES);
	  glColor3f(1,0,0);
	  glVertex3f(0,0,0);
	  glVertex3f(40,0,0);

	  glColor3f(0,1,0);
	  glVertex3f(0,0,0);
	  glVertex3f(0,40,0);

	  glColor3f(0,0,1);
	  glVertex3f(0,0,0);
	  glVertex3f(0,0,40);
     glEnd();

     // draw a sphere 
     if ( !datafile ) 
     {
       glColor3f(0.9, 0.7, 0.1);
       glutWireSphere(30, 10, 10);
     }

   glPopMatrix();

   // draw a box
   glColor3f(1.0, 1.0, 1.0);
   glPushMatrix();

   // box behind the screen when loaded
   glTranslatef( - (Xminmax[0] + Xminmax[1])/2.0, 
		   - (Yminmax[0] + Yminmax[1])/2.0, 
		   - Zminmax[1] );

   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

   glBegin(GL_POLYGON);
	   glVertex3f(Xminmax[1], Yminmax[1], Zminmax[0]); // 1
	   glVertex3f(Xminmax[1], Yminmax[0], Zminmax[0]); // 2
	   glVertex3f(Xminmax[0], Yminmax[0], Zminmax[0]); // 3 
	   glVertex3f(Xminmax[0], Yminmax[1], Zminmax[0]); // 4
   glEnd();

   glBegin(GL_POLYGON);
	   glVertex3f(Xminmax[1], Yminmax[1], Zminmax[0]); // 1
	   glVertex3f(Xminmax[1], Yminmax[1], Zminmax[1]); // 8
	   glVertex3f(Xminmax[1], Yminmax[0], Zminmax[1]); // 7 
	   glVertex3f(Xminmax[1], Yminmax[0], Zminmax[0]); // 2
   glEnd();

   glBegin(GL_POLYGON);
	   glVertex3f(Xminmax[1], Yminmax[1], Zminmax[0]); // 1
	   glVertex3f(Xminmax[0], Yminmax[1], Zminmax[0]); // 4
	   glVertex3f(Xminmax[0], Yminmax[1], Zminmax[1]); // 5
	   glVertex3f(Xminmax[1], Yminmax[1], Zminmax[1]); // 8
   glEnd();

   glBegin(GL_POLYGON);
	   glVertex3f(Xminmax[0], Yminmax[1], Zminmax[0]); // 4
	   glVertex3f(Xminmax[0], Yminmax[0], Zminmax[0]); // 3 
	   glVertex3f(Xminmax[0], Yminmax[0], Zminmax[1]); // 6
	   glVertex3f(Xminmax[0], Yminmax[1], Zminmax[1]); // 5
   glEnd();

   glBegin(GL_POLYGON);
	   glVertex3f(Xminmax[1], Yminmax[0], Zminmax[0]); // 2
	   glVertex3f(Xminmax[1], Yminmax[0], Zminmax[1]); // 7 
	   glVertex3f(Xminmax[0], Yminmax[0], Zminmax[1]); // 6
	   glVertex3f(Xminmax[0], Yminmax[0], Zminmax[0]); // 3 
   glEnd();

   glBegin(GL_POLYGON);
	   glVertex3f(Xminmax[0], Yminmax[1], Zminmax[1]); // 5
	   glVertex3f(Xminmax[0], Yminmax[0], Zminmax[1]); // 6
	   glVertex3f(Xminmax[1], Yminmax[0], Zminmax[1]); // 7 
	   glVertex3f(Xminmax[1], Yminmax[1], Zminmax[1]); // 8
   glEnd();

   glPopMatrix();

   // draw the brain dataset
   if ( datafile ) 
   {
     glPushMatrix();

     // model behind the screen when loaded
     glTranslatef( 0., 0.,  -(Zminmax[0]+Zminmax[1]) /2.);
     glRotatef(rot_x, 1,0,0);
     glRotatef(rot_y, 0,1,0);
     glRotatef(rot_z, 0,0,1);

     tubes.draw();

     glPopMatrix();
   }

   glutSwapBuffers();
}

/*------------------------------------------------------------------------
 * Display call back
 *
 *  Input: none
 *  Output:
 *   it will draw a box (which is the same size as the brain) and a yellow
 *   sphere, if no datafile is specified (when running ./stereoTest -stereo)
 *
 *   Otherwise, it will draw a box and a brain. 
 *   (when running ./stereoTest -stereo -f ./test.data)
 *
 *   The model is set to be behind the screen when it is first loaded.
 */
static void draw(void)
{
    XYZ r;
    double ratio, radians, wd2, ndfl;
    double left, right, top, bottom;

    /* Clip to avoid extreme stereo */
    double near = camera.focallength / 5;
    double far = (Yminmax[1]-Yminmax[2]) * 40.0;

    /* Misc stuff */
    ratio = camera.screenwidth / (double) camera.screenheight;
    radians = DTOR * camera.aperture / 2.;
    wd2 = near * tan(radians);
    ndfl = near / camera.focallength;

    /* Clear the buffers */
    glDrawBuffer(GL_BACK_LEFT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawBuffer(GL_BACK_RIGHT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /* Derive the two eye positions */
    CROSSPROD(camera.vd, camera.vu, r);
    Normalise(&r);
    r.x *= camera.eyesep / 2.0;
    r.y *= camera.eyesep / 2.0;
    r.z *= camera.eyesep / 2.0;

    // left eye
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        left = -ratio * wd2 - 0.5 * camera.eyesep * ndfl;
	right = ratio * wd2 - 0.5 * camera.eyesep * ndfl;
	top = wd2;
	bottom = -wd2;
	glFrustum(left, right, bottom, top, near*2, far);

	glMatrixMode(GL_MODELVIEW);
	glDrawBuffer(GL_BACK_RIGHT);
	glLoadIdentity();
	gluLookAt(camera.vp.x + r.x, camera.vp.y + r.y, camera.vp.z + r.z,
		camera.vp.x + r.x + camera.vd.x,
		camera.vp.y + r.y + camera.vd.y,
		camera.vp.z + r.z + camera.vd.z,
		camera.vu.x, camera.vu.y, camera.vu.z);
        do_draw();

    // right eye
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        left = -ratio * wd2 + 0.5 * camera.eyesep * ndfl;
	right = ratio * wd2 + 0.5 * camera.eyesep * ndfl;
	top = wd2;
	bottom = -wd2;
	glFrustum(left, right, bottom, top, near*2, far);

	glMatrixMode(GL_MODELVIEW);
	glDrawBuffer(GL_BACK_LEFT);
	glLoadIdentity();
	gluLookAt(camera.vp.x - r.x, camera.vp.y - r.y, camera.vp.z - r.z,
				camera.vp.x - r.x + camera.vd.x,
				camera.vp.y - r.y + camera.vd.y,
				camera.vp.z - r.z + camera.vd.z,
				camera.vu.x, camera.vu.y, camera.vu.z);
	do_draw();

    //glutSwapBuffers();
}

/*------------------------------------------------------------------------
 * resize window 
 */
static void reshape(int width, int height)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    camera.screenwidth = width;
    camera.screenheight = height;
}

/*------------------------------------------------------------------------
 * OpenGL main routine
 */
int main(int argc, char *argv[])
{
   int i;

   for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-f") == 0) {
		  if ( i+1 >= argc ) {
			  printf("-f: missing data file.\n");
			  return 1;
		  }
		  datafile = GL_TRUE;
		  break;
      }
      else {
		 printf("Warrning: unknown parameter: %s\n", argv[i]);
	  }
   }

   if ( datafile ) {
	   if ( 0 != tubes.loadGeometry( argv[argc-1] ) ) {
		   return 1;
	   }
	   tubes.getExtent( Xminmax, Yminmax, Zminmax );
   }

   glutInit (&argc, argv);
   glutInitWindowPosition ( 0, 0 );
   glutInitWindowSize ( 1920, 1080 );
   glutInitDisplayMode	( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO );
   glutCreateWindow ("brain stereoscopic rendering for visweek");
   initCamera();
   glutReshapeFunc(reshape);
   glutDisplayFunc(draw);
   glutKeyboardFunc(keyFunc);
   glutMainLoop();

   return 0;
}

