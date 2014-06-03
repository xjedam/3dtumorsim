#include <inttypes.h>
#include <GL/glut.h>
#include <stdio.h>
#include <time.h>

#define WINDOW_NAME       "3d tumor simulator"
#define MODEL_SIZE_X      30
#define MODEL_SIZE_Y      30
#define MODEL_SIZE_Z      30
#define WINDOW_X          1000
#define WINDOW_Y          1000
#define SITE_SIZE         7.0f
#define ITER_DELAY        0.2
#define DEBUG(a)          a

#if defined(__CYGWIN__) || defined(_WIN64) || defined(_WIN32)
  #include <windows.h>
  #define SLEEP_FUNC Sleep
	#define SLEEP_MULTIPLIER 1000
#else
	#include <unistd.h>
  #define SLEEP_FUNC usleep
	#define SLEEP_MULTIPLIER 1000000
#endif

extern int numCells;
  
void processNormalKeys(unsigned char key, int x, int y);
void processSpecialKeys(int key, int x, int y);
void onReshape(int w, int h);
void display();
