#include "3dtumorsim.h"
#include "lattice.h"
#include "potts.h"

// main sub-cell latice
int64_t ***lattice;

// cell array
cell_info_t *cells;

// rotation
float xrot = 0.0, yrot = 0.0;

int xClick, yClick, width, height;
int lmbDown = 0;
int isPause = 1;
int numCells = 0;
clock_t start, stop;
int iterationNumber = 0;

// initialize OpenGL perspective
void initGL() {
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f);                // Set background color to black and opaque
   glClearDepth(1.0f);                                  // Set background depth to farthest
   glEnable(GL_DEPTH_TEST);                             // Enable depth testing for z-culling
   glDepthFunc(GL_LEQUAL);                              // Set the type of depth-test
   glShadeModel(GL_SMOOTH);                             // Enable smooth shading
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);   // Nice perspective corrections
}

// Clears the current window and draws the iteration.
void display() {
  if(!isPause) {
    stop = clock();
    if(((float)(stop - start))/CLOCKS_PER_SEC < ITER_DELAY && !isPause){
      SLEEP_FUNC((ITER_DELAY - (((float)(stop - start))/CLOCKS_PER_SEC))*SLEEP_MULTIPLIER);
    }

    start = clock();
    
    iterationNumber++;
    printf("Iteration: %i, number of cells: %i\n", iterationNumber, numCells);
    calculateNextStep(lattice, numCells, cells);
  }

  glClearDepth(1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);             // To operate on model-view matrix
  drawCells(xrot, yrot, lattice, cells);
  glutSwapBuffers();
}

// Handles keypress
void processNormalKeys(unsigned char key, int x, int y) {
	if (key == 27) {
		exit(0);
  } else if(key == 'p') {
    isPause = isPause == 0;
  } 
}

// Handles mouse press
void mousePressCallback(int button, int state, int x, int y) {
  if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    lmbDown = 1;
    xClick = x;
    yClick = y;
  } else if(button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
    lmbDown = 0;
  }
}

// Handles mouse drag
void mouseActiveMoveCallback(int x, int y) {
  if(lmbDown) {
    xrot += xClick - x;
    yrot += yClick - y;  
    xClick = x;
    yClick = y;
    onReshape(width, height);
  }
}

// recalculates view on window resize
void onReshape(int w, int h) {
  width = w;
  height = h;
	glMatrixMode(GL_PROJECTION); // projection matrix is active
	glLoadIdentity(); // reset the projection

  // set rotation and camera position
  glRotatef(yrot, 1.0f, 0.0f, 0.0f);
  glRotatef(xrot, 0.0f, 1.0f, 0.0f);
  glFrustum(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5);

	glMatrixMode(GL_MODELVIEW); // return to modelview mode
  glViewport(SITE_SIZE/2, SITE_SIZE/2, w, h);
}


int main(int argc, char** argv) {

  srand(time(NULL));

  // init lattioce and cells array
  lattice = initLattice();
  cells = initCells();

  // Use a single buffered window in RGB mode (as opposed to a double-buffered
  // window or color-index mode).
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  start = clock();
  
  // Position window at (80,80) and give it a title.
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
  glutReshapeFunc(onReshape);       
  glutMotionFunc(mouseActiveMoveCallback);
  glutMouseFunc(mousePressCallback);
  
  glutKeyboardFunc(processNormalKeys);
  
  if(argc > 1) {
    FILE *fp = fopen(argv[1], "r");
    numCells = loadModel(fp, lattice, cells);
  }

  // Tell GLUT to start reading and processing events.  This function
  // never returns; the program only exits when the user closes the main
  // window or kills the process.
  glutMainLoop();

  return 0;
}