#define	MEMBRANE_RESISTANCE 20
#define VOLUME_RESISTANCE 20
#define ENERGY_TRESHOLD 20
#define TEMPERATURE 36.7

typedef struct POTTS_T {
  int membraneArea;
  int volume;
  int type;
  float temperature;
} cell_info_t;

cell_info_t *initCells();
int loadModel(FILE *fp, int64_t ***lattice, cell_info_t *cells);
void saveModel(FILE *fp, int64_t ***lattice, cell_info_t *cells);
int getMembraneChangeInsert(int x, int y, int z, int64_t ***lattice, int sigma);
int getMembraneChangeRemove(int x, int y, int z, int64_t ***lattice, int sigma);
void calculateNextStep(int64_t ***lattice, int cellCount, cell_info_t *cells);
int checkAndMoveToSite(int x1, int y1, int z1, int x2, int y2, int z2, int64_t ***lattice, cell_info_t *cells);