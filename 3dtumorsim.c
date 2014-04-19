#include "3dtumorsim.h"
#include "lattice.h"
#include <time.h>

int64_t ***lattice, ***buff;
int isPause = 0;
clock_t start, stop;

// Clears the current window and draws a triangle.
void display() {
  int i, j, k;
  if(!isPause) {
    stop = clock();
    if(((float)(stop - start))/CLOCKS_PER_SEC < ITER_DELAY && !isPause){
      //printf("%f\n", ITER_DELAY - (((float)(stop - start))/CLOCKS_PER_SEC));
      SLEEP_FUNC((ITER_DELAY - (((float)(stop - start))/CLOCKS_PER_SEC))*SLEEP_MULTIPLIER);
    }
    start = clock();
    // Set every pixel in the frame buffer to the current clear color.
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on model-view matrix
    
    for(i = 0; i < MODEL_SIZE_X; i++) {
      for(j = 0; j < MODEL_SIZE_Y; j++) {
        for(k = 0; k < MODEL_SIZE_Z; k++) {
          if(lattice[i][j][k] != 0) {
            drawLatticeSite(i, j, k, lattice[i][j][k]);
          }
        }
      }
    }
    
    // Drawing is done by specifying a sequence of vertices.  The way these
    // vertices are connected (or not connected) depends on the argument to
    // glBegin.  GL_POLYGON constructs a filled polygon.
    // glBegin(GL_POLYGON);
      // glColor3f(1, 0, 0); glVertex3f(-0.6, -0.75, 0.5);
      // glColor3f(0, 1, 0); glVertex3f(0.6, -0.75, 0);
      // glColor3f(0, 0, 1); glVertex3f(0, 0.75, 0);
    // glEnd();

    glutSwapBuffers();
    // Flush drawing command buffer to make drawing happen as soon as possible.
    // glFlush();
  }
}

void processNormalKeys(unsigned char key, int x, int y) {
	if (key == 27) {
		exit(0);
  } else if(key == 'p') {
    isPause = isPause == 0;
  } else if(key == 'p') {
    isPause = isPause == 0;
  } else if(key == 'p') {
    isPause = isPause == 0;
  }
}

void processSpecialKeys(int key, int x, int y) {

	// switch(key) {
		// case GLUT_KEY_F1 :
				// red = 1.0;
				// green = 0.0;
				// blue = 0.0; break;
		// case GLUT_KEY_F2 :
				// red = 0.0;
				// green = 1.0;
				// blue = 0.0; break;
		// case GLUT_KEY_F3 :
				// red = 0.0;
				// green = 0.0;
				// blue = 1.0; break;
	// }
}

// Initializes GLUT, the display mode, and main window; registers callbacks;
// enters the main event loop.
int main(int argc, char** argv) {

  lattice = initLattice();
  buff = initLattice();
  // Use a single buffered window in RGB mode (as opposed to a double-buffered
  // window or color-index mode).
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  start = clock();
  
  // Position window at (80,80)-(480,380) and give it a title.
  glutInitWindowPosition(80, 80);
  glutInitWindowSize(WINDOW_X, WINDOW_Y);
  glutCreateWindow(WINDOW_NAME);
  
  // Enable alpha blending
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable( GL_BLEND );
  // Tell GLUT that whenever the main window needs to be repainted that it
  // should call the function display().
  glutDisplayFunc(display);
  glutIdleFunc(display);
  // glutReshapeFunc(reshape);       // Register callback handler for window re-size event
  
  glutKeyboardFunc(processNormalKeys);
	glutSpecialFunc(processSpecialKeys);
  
  int i = 0;
  for(i = 0; i < MODEL_SIZE_X; i++) {
    lattice[i][i][i] = SET_TYPE(0, 1);
  }

  // Tell GLUT to start reading and processing events.  This function
  // never returns; the program only exits when the user closes the main
  // window or kills the process.
  glutMainLoop();
}