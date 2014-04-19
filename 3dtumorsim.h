#include <inttypes.h>
#include <GL/glut.h>

#define WINDOW_NAME       "3d tumor simulator"
#define MODEL_SIZE_X      50
#define MODEL_SIZE_Y      50
#define MODEL_SIZE_Z      50
#define WINDOW_X          500
#define WINDOW_Y          500
#define SITE_SIZE         7
#define ITER_DELAY        0.4

#if defined(__CYGWIN__) || defined(_WIN64) || defined(_WIN32)
  #include <windows.h>
  #define SLEEP_FUNC Sleep
	#define SLEEP_MULTIPLIER 1000
#else
	#include <unistd.h>
  #define SLEEP_FUNC usleep
#define SLEEP_MULTIPLIER 1000000
#endif

void processNormalKeys(unsigned char key, int x, int y);
void processSpecialKeys(int key, int x, int y);
void display();