#define	MEMBRANE_RESISTANCE 10
#define VOLUME_RESISTANCE 10
#define ENERGY_TRESHOLD 100
#define TEMPERATURE 10
#define CELL_OVERSIZE_SPLIT_TRESHOLD 0.30			

typedef struct SITE {
	int x;
	int y;
	int z;
	int64_t *val;
} site_t;

typedef struct POTTS_T {
  int membraneArea;
  int volume;
  int type;
  float temperature;
  site_t *subcells;
} cell_info_t;

cell_info_t *initCells();
int loadModel(FILE *fp, int64_t ***lattice, cell_info_t *cells);
void saveModel(FILE *fp, int64_t ***lattice, cell_info_t *cells);
int getMembraneChangeInsert(int x, int y, int z, int64_t ***lattice, int sigma);
int getMembraneChangeRemove(int x, int y, int z, int64_t ***lattice, int sigma);
int getAdhesionChangeInsert(int x, int y, int z, int64_t ***lattice, int sigma);
int getAdhesionChangeRemove(int x, int y, int z, int64_t ***lattice, int sigma);
void calculateNextStep(int64_t ***lattice, int cellCount, cell_info_t *cells);
int checkAndMoveToSite(int x1, int y1, int z1, int x2, int y2, int z2, int64_t ***lattice, cell_info_t *cells);
void splitCell(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigma);
long long getNewEnergyLatticeChange(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigmaTo, int sigmaFrom);
void changeSiteOwner(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigmaTo, int sigmaFrom);