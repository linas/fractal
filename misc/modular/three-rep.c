
/* 
 * Draw orbifold for the 3D representation of the
 * modular group.
 *
 * November 2004
 * Copyright (c) 2004 Linas Vepstas <linas@linas.org>
 */


/* required include files */
#include <stdlib.h>

#include <GL/gl.h>
#include <GL/glut.h>

extern void InitStuff(void);
extern void DrawStuff(void);

float lastx=0;
float lasty=0;

/* get notified of mouse motions */
static void 
MouseMotion (int x, int y)
{
   lastx = x;
   lasty = y;
   glutPostRedisplay ();
}

static void 
InitMouse (void) 
{
   lastx = 121.0;
   lasty = 121.0;
}

/* draw the helix shape */
void 
DrawSpin (void) 
{

   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glColor3f (0.6, 0.8, 0.3);

   /* set up some matrices so that the object spins with the mouse */
   glPushMatrix ();
   glTranslatef (0.0, 0.0, -80.0);
   glRotatef (lastx, 0.0, 1.0, 0.0);
   glRotatef (lasty, 1.0, 0.0, 0.0);

   /* Phew. FINALLY, stuff  -- */
	DrawStuff();

   glPopMatrix ();

   glutSwapBuffers ();
}

/* ------------------------- end of file ----------------- */


/* set up a light */
GLfloat lightOnePosition[] = {40.0, 40, 100.0, 0.0};
GLfloat lightOneColor[] = {0.99, 0.99, 0.99, 1.0}; 

GLfloat lightTwoPosition[] = {-40.0, 40, 100.0, 0.0};
GLfloat lightTwoColor[] = {0.99, 0.99, 0.99, 1.0}; 

int
main (int argc, char * argv[]) 
{
   /* initialize glut */
   glutInit (&argc, argv);
   glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   glutCreateWindow ("join styles");
   glutDisplayFunc (DrawSpin);
   glutMotionFunc (MouseMotion);

   /* initialize GL */
   glClearDepth (1.0);
   glEnable (GL_DEPTH_TEST);
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glShadeModel (GL_SMOOTH);

   glMatrixMode (GL_PROJECTION);
   /* roughly, measured in centimeters */
   glFrustum (-9.0, 9.0, -9.0, 9.0, 50.0, 150.0);
   glMatrixMode(GL_MODELVIEW);

#if 0
   /* initialize lighting */
   glLightfv (GL_LIGHT0, GL_POSITION, lightOnePosition);
   glLightfv (GL_LIGHT0, GL_DIFFUSE, lightOneColor);
   glEnable (GL_LIGHT0);
   glLightfv (GL_LIGHT1, GL_POSITION, lightTwoPosition);
   glLightfv (GL_LIGHT1, GL_DIFFUSE, lightTwoColor);
   glEnable (GL_LIGHT1);
   glEnable (GL_LIGHTING);
   glEnable (GL_NORMALIZE);
   glColorMaterial (GL_FRONT_AND_BACK, GL_DIFFUSE);
   glEnable (GL_COLOR_MATERIAL);
#endif

   InitMouse ();
   InitStuff ();

   glutMainLoop ();
   return 0;             /* ANSI C requires main to return int. */
}
/* ------------------ end of file -------------------- */
