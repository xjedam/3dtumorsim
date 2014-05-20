#include "3dtumorsim.h"
#include "lattice.h"
#include "potts.h"

int64_t ***lattice, ***buff;
cell_info_t *cells;
float xrot = 0.0, yrot = 0.0;
int xClick, yClick, width, height;
int lmbDown = 0;
int isPause = 0;
int numCells = 0;
clock_t start, stop;

void initGL() {
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
   glClearDepth(1.0f);                   // Set background depth to farthest
   glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
   glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
   glShadeModel(GL_SMOOTH);   // Enable smooth shading
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections
}

// Clears the current window and draws a triangle.
void display() {
  int i, j, k;
  DEBUG(printf("rot: %f, %f\n", xrot, yrot);)
  if(!isPause) {
    stop = clock();
    if(((float)(stop - start))/CLOCKS_PER_SEC < ITER_DELAY && !isPause){
      SLEEP_FUNC((ITER_DELAY - (((float)(stop - start))/CLOCKS_PER_SEC))*SLEEP_MULTIPLIER);
    }
    start = clock();
    // Set every pixel in the frame buffer to the current clear color.
    glClearDepth(1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on model-view matrix
    
    // for(i = 0; i < MODEL_SIZE_X; i++) {
    //   for(j = 0; j < MODEL_SIZE_Y; j++) {
    //     for(k = 0; k < MODEL_SIZE_Z; k++) {
    //       if(lattice[i][j][k] != 0) {
    //         drawLatticeSite(i, j, k, lattice[i][j][k], xrot, yrot, lattice);
    //       }
    //     }
    //   }
    // }

    drawCells(xrot, yrot, lattice, cells);

    calculateNextStep(lattice, numCells, cells);
    
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

void mousePressCallback(int button, int state, int x, int y) {
  if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    lmbDown = 1;
    xClick = x;
    yClick = y;
  } else if(button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
    lmbDown = 0;
  }
}

void mouseActiveMoveCallback(int x, int y) {
  if(lmbDown) {
    xrot += xClick - x;
    yrot += yClick - y;  
    xClick = x;
    yClick = y;
    onReshape(width, height);
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

void onReshape(int w, int h) {
  width = w;
  height = h;
	glMatrixMode(GL_PROJECTION); // projection matrix is active
	glLoadIdentity(); // reset the projection
  glRotatef(yrot, 1.0f, 0.0f, 0.0f);
  glRotatef(xrot, 0.0f, 1.0f, 0.0f);
  glFrustum(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5);
	// gluPerspective(45.0, ratio, 0.1, 100.0); // perspective transformation
	glMatrixMode(GL_MODELVIEW); // return to modelview mode
  glViewport(SITE_SIZE/2, SITE_SIZE/2, w, h);
}

// Initializes GLUT, the display mode, and main window; registers callbacks;
// enters the main event loop.
int main(int argc, char** argv) {

  srand(time(NULL));
  lattice = initLattice();
  buff = initLattice();
  cells = initCells();
  // Use a single buffered window in RGB mode (as opposed to a double-buffered
  // window or color-index mode).
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  start = clock();
  
  // Position window at (80,80)-(480,380) and give it a title.
  glutInitWindowPosition(80, 80);
  glutInitWindowSize(WINDOW_X, WINDOW_Y);
  glutCreateWindow(WINDOW_NAME);
  
  // Enable alpha blending
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable( GL_BLEND );
  
  initGL();
  
  // Tell GLUT that whenever the main window needs to be repainted that it
  // should call the function display().
  glutDisplayFunc(display);
  glutIdleFunc(display);
  glutReshapeFunc(onReshape);       // Register callback handler for window re-size event
  glutMotionFunc(mouseActiveMoveCallback);
  glutMouseFunc(mousePressCallback);
  
  glutKeyboardFunc(processNormalKeys);
	glutSpecialFunc(processSpecialKeys);
  
  if(argc > 1) {
    FILE *fp = fopen(argv[1], "r");
    numCells = loadModel(fp, lattice, cells);
  }
  // int i = 0;
  // for(i = 0; i < MODEL_SIZE_X; i++) {
  //   lattice[i][i][i] = SET_TYPE(1, (VASCULAR + i) % 2 + 1 );
  //   printf("placing at [%i, %i, %i] %i\n", i, i, i, (VASCULAR + i) % 2 + 1);
  // }

  // Tell GLUT to start reading and processing events.  This function
  // never returns; the program only exits when the user closes the main
  // window or kills the process.
  glutMainLoop();

  return 0;
}