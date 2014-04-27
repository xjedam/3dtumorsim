typedef struct POTTS_T {
  int membraneArea;
  int volume;
  int type;
  float temperature;
} cell_info_t;

cell_info_t *initCells();
int loadModel(FILE *fp, int64_t ***lattice, cell_info_t *cells);
void saveModel(FILE *fp, int64_t ***lattice, cell_info_t *cells);
int getMembraneChange(int x, int y, int z, int64_t ***lattice);