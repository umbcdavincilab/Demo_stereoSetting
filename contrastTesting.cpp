/*
 *
 *  Filename:  contrastTesting.cpp
 *
 *  Description:
 *    This code is used to show there patterns: 
 *           - full black
 *           - full white
 *           - ANSI 4x4 block patterns
 *
 *    The code is written on a Linux platform.
 *
 *  Usage:
 *    ./contrastTesting
 *            use default screen size 1920x1028
 *    ./contrastTesting [width] [height]
 *            use the user-defined screen size
 *  
 *  Creation:  May 25 2014
 *  Author: Jian Chen 
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>
#include "GLoader.h"

int next_image = 0;  // for shuffling the three images

int width = 1920, height = 1028; // default screen size

/*------------------------------------------------------------------------
 *  init
 *    Input:
 *    Output:
 *    Description:
 *       init env
 *  
 */
static void init()
{
  glClearColor(0,0,0,0);
  glShadeModel(GL_FLAT);
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
		case 'n':
		case 'N':
		    next_image=next_image++;
			if(next_image==3) next_image=0;
			break;
		default:
			return;
	}
	glutPostRedisplay();
}


/*------------------------------------------------------------------------
 * draw the checkboard pattern
 *
 *  Input: none
 *  Output:
 */
void drawChecker(int w_size, int h_size)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < 4; i++) {
	  for( j=0; j<4; j++)
	  {
        if((i + j)%2 == 0) // if i + j is even
            glColor3f( 1, 1, 1);
        else
            glColor3f( 0, 0, 0);
        glRecti(i*w_size, j*h_size, (i+1)*w_size, (j+1)*h_size);    // draw the rectangle
      }; // end for(j)
    }; // end for(i)
}


/*------------------------------------------------------------------------
 * Display call back
 *
 *  Input: none
 *  Output:
 */
static void draw(void)
{
  glClear (GL_COLOR_BUFFER_BIT);
  glColor3f(1,1,1);

  glDisable(GL_TEXTURE_2D);

  if(next_image==0)  // all white
  {
    glColor3f(1,1,1);
    glRecti(0, 0, width, height);    // draw the rectangle
  }
  else if(next_image==1) // all black
  {
    glColor3f(0,0,0);
    glRecti(0, 0, width, height);    // draw the rectangle
  }
  else if(next_image==2) // ansi 4x4
  {
    drawChecker(width/4, height/4);
  }
  glFlush();
}

/*------------------------------------------------------------------------
 * resize window 
 */
static void reshape(int width, int height)
{
   glViewport (0, 0, (GLsizei) width, (GLsizei) height);
   glMatrixMode (GL_PROJECTION);
   glLoadIdentity ();
   gluOrtho2D (0.0, (GLdouble) width, 0.0, (GLdouble) height);
}

/*------------------------------------------------------------------------
 * OpenGL main routine
 */
int main(int argc, char *argv[])
{
   int i;
   
   if(argc!=3) 
   {
     printf("Use default screen size: 1920x1028. \n");
   }
   else if(argc==3)
   {
     width = atoi(argv[1]); 
	 height = atoi(argv[2]);
   }

   glutInit (&argc, argv);
   glutInitWindowPosition ( 0, 0 );
   glutInitWindowSize ( width, height);
   glutInitDisplayMode	( GLUT_SINGLE| GLUT_RGB );
   glutCreateWindow ("screen contrast measurement patterns for visweek");
   init();
   glutReshapeFunc(reshape);
   glutDisplayFunc(draw);
   glutKeyboardFunc(keyFunc);
   glutMainLoop();

   return 0;
}

