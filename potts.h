#define	MEMBRANE_RESISTANCE 20
#define VOLUME_RESISTANCE 20

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
void calculateNextStep(int64_t ***lattice, int cellCount);