#include "3dtumorsim.h"
#include "lattice.h"
#include "potts.h"
#include <math.h>

// defining bond energy
int bondEnergy[6][6] = {
  {0, 12, 10, 15, 10, 0},       // MEDIUM
  {12, 160, 50, 30, 30, 0},      // VASCULAR
  {10, 50, 90, 8, 8, 0},        // TUMOR_NORM
  {15, 30, 8, 40, 8, 0},        // TUMOR_NECROSIS
  {10, 30, 8, 8, 8, 0},         // TUMOR_STEM
  {0, 0, 0, 0, 0, 0},
};

// target membrane area definition
int targetMembrane[6] = {
  30,                           // MEDIUM
  55,                           // VASCULAR
  50,                           // TUMOR_NORM
  50,                           // TUMOR_NECROSIS
  50,                           // TUMOR_STEM
  0
};

// target volume definition
int targetVolume[6] = {
  15,                           // MEDIUM
  25,                           // VASCULAR
  15,                           // TUMOR_NORM
  15,                           // TUMOR_NECROSIS
  15,                           // TUMOR_STEM
  0
};

energy_t globalEnergy = {0, 0, 0, 0};

cell_info_t *initCells() {
	int i, size = MODEL_SIZE_X * MODEL_SIZE_Y * MODEL_SIZE_Z / 8;
	cell_info_t *cells = (cell_info_t *)malloc(size * sizeof(cell_info_t));
	for(i = 0; i < size; i++) {
		cells[i].temperature = 36.7;
		cells[i].type = MEDIUM;
		cells[i].volume = 15;
		cells[i].membraneArea = 30;
	}
	return cells;
}

void saveModel(FILE *fp, int64_t ***lattice, cell_info_t *cells) {
  int i, j, k;
  for(i = 0; i < MODEL_SIZE_X; i++) {
    for(j = 0; j < MODEL_SIZE_Y; j++) {
      for(k = 0; k < MODEL_SIZE_Z; k++) {
        int site = lattice[i][j][k];
        if(TYPE(site) != MEDIUM) {
          fprintf(fp, "%i %i %i: %i %i %f\n", i, j, k, TYPE(site), SIGMA(site), cells[SIGMA(site)].temperature);
        }
      }
    }
  }

  fclose(fp);
}

// returns ammount of cells loaded
int loadModel(FILE *fp, int64_t ***lattice, cell_info_t *cells) {
	int x, y, z, sigma, type, i, max = 0, j, k;
  float temperature;
	while(fscanf(fp, "%i%i%i:%i%i%f", &x, &y, &z, &type, &sigma, &temperature) != EOF) {
    if(cells[sigma].type != type) {
      cells[sigma].subcells = (site_t *)malloc(targetVolume[type] * 2 * sizeof(site_t));
      cells[sigma].type = type;
      cells[sigma].temperature = temperature;
      cells[sigma].volume = 0;
      cells[sigma].membraneArea = 0;
    }
    cells[sigma].membraneArea += getMembraneChangeInsert(x, y, z, lattice, sigma);
		lattice[x][y][z] = SET_TYPE(SET_SIGMA(0,sigma), type);
    cells[sigma].subcells[cells[sigma].volume].val = &lattice[x][y][z];
    cells[sigma].subcells[cells[sigma].volume].x = x;
    cells[sigma].subcells[cells[sigma].volume].y = y;
    cells[sigma].subcells[cells[sigma].volume].z = z;
		cells[sigma].volume += 1;
    if(sigma > max) {
      max = sigma;
    }
	}

  // calculate global energy: volume, membrane
  int membranePart = 0;
  int volumePart = 0;
  for(i = 1; i <= max; i++) {
    membranePart += pow(cells[i].membraneArea - targetMembrane[cells[i].type], 2);
    volumePart += pow(cells[i].volume - targetVolume[cells[i].type], 2);
  }
  globalEnergy.membrane = MEMBRANE_RESISTANCE * membranePart;
  globalEnergy.volume = VOLUME_RESISTANCE * volumePart;

  // add global energy: adhesion part
  for(i = 0; i < MODEL_SIZE_X; i++) {
    for(j = 0; j < MODEL_SIZE_Y; j++) {
      for(k = 0; k < MODEL_SIZE_Z; k++) {
        if(i + 1 < MODEL_SIZE_X && SIGMA(lattice[i][j][k]) != SIGMA(lattice[i + 1][j][k])) {
          globalEnergy.adhesion += bondEnergy[SIGMA(lattice[i][j][k])][SIGMA(lattice[i + 1][j][k])];
        }
        if(j + 1 < MODEL_SIZE_Y && SIGMA(lattice[i][j][k]) != SIGMA(lattice[i][j + 1][k])) {
          globalEnergy.adhesion += bondEnergy[SIGMA(lattice[i][j][k])][SIGMA(lattice[i][j + 1][k])];
        }
        if(k + 1 < MODEL_SIZE_Z && SIGMA(lattice[i][j][k]) != SIGMA(lattice[i][j][k + 1])) {
          globalEnergy.adhesion += bondEnergy[SIGMA(lattice[i][j][k])][SIGMA(lattice[i][j][k + 1])];
        }
      }
    }
  }
  globalEnergy.total = globalEnergy.membrane + globalEnergy.volume + globalEnergy.adhesion;

	fclose(fp);
  DEBUG(puts("loaded model");)
  return max;
}

// gets the change of the membrane area when x, y, z was inserted
int getMembraneChangeInsert(int x, int y, int z, int64_t ***lattice, int sigma) {
  int change = -6;

  if(x == 0 || SIGMA(lattice[x - 1][y][z]) != sigma) {
    change += 2;
  }
  if(x == MODEL_SIZE_X - 1 || SIGMA(lattice[x + 1][y][z]) != sigma) {
    change += 2;
  }
  if(y == 0 || SIGMA(lattice[x][y - 1][z]) != sigma) {
    change += 2;
  }
  if(y == MODEL_SIZE_Y - 1 || SIGMA(lattice[x][y + 1][z]) != sigma) {
    change += 2;
  }
  if(z == 0 || SIGMA(lattice[x][y][z - 1]) != sigma) {
    change += 2;
  }
  if(z == MODEL_SIZE_Z - 1 || SIGMA(lattice[x][y][z + 1]) != sigma) {
    change += 2;
  }

  return change;
}

// gets the change of the membrane area when x, y, z was removed
int getMembraneChangeRemove(int x, int y, int z, int64_t ***lattice, int sigma) {
  int change = -6;

  if(x > 0 && SIGMA(lattice[x - 1][y][z]) == sigma) {
    change += 2;
  }
  if(x < MODEL_SIZE_X - 1 && SIGMA(lattice[x + 1][y][z]) == sigma) {
    change += 2;
  }
  if(y > 0 && SIGMA(lattice[x][y - 1][z]) == sigma) {
    change += 2;
  }
  if(y < MODEL_SIZE_Y - 1 && SIGMA(lattice[x][y + 1][z]) == sigma) {
    change += 2;
  }
  if(z > 0 && SIGMA(lattice[x][y][z - 1]) == sigma) {
    change += 2;
  }
  if(z < MODEL_SIZE_Z - 1 && SIGMA(lattice[x][y][z + 1]) == sigma) {
    change += 2;
  }

  return change;
}

int getAdhesionChangeRemove(int x, int y, int z, int64_t ***lattice, int sigma) {
  int change = 0, type = TYPE(lattice[x][y][z]);

  if(x > 0 && SIGMA(lattice[x - 1][y][z]) != sigma) {
    change -= bondEnergy[TYPE(lattice[x - 1][y][z])][type];
  }
  if(x < MODEL_SIZE_X - 1 && SIGMA(lattice[x + 1][y][z]) != sigma) {
    change -= bondEnergy[TYPE(lattice[x + 1][y][z])][type];
  }
  if(y > 0 && SIGMA(lattice[x][y - 1][z]) != sigma) {
    change -= bondEnergy[TYPE(lattice[x][y - 1][z])][type];
  }
  if(y < MODEL_SIZE_Y - 1 && SIGMA(lattice[x][y + 1][z]) != sigma) {
    change -= bondEnergy[TYPE(lattice[x][y + 1][z])][type];
  }
  if(z > 0 && SIGMA(lattice[x][y][z - 1]) != sigma) {
    change -= bondEnergy[TYPE(lattice[x][y][z - 1])][type];
  }
  if(z < MODEL_SIZE_Z - 1 && SIGMA(lattice[x][y][z + 1]) != sigma) {
    change -= bondEnergy[TYPE(lattice[x][y][z + 1])][type];
  }

  DEBUG(printf("\t\t\tAdhesion remove from [%i, %i, %i] sig %i, %i\n", x, y, z, sigma, change);)
  return change;
}

int getAdhesionChangeInsert(int x, int y, int z, int64_t ***lattice, int sigma) {
  int change = 0, type = TYPE(lattice[x][y][z]);

  if(x > 0 && SIGMA(lattice[x - 1][y][z]) != sigma) {
    change += bondEnergy[TYPE(lattice[x - 1][y][z])][type];
  }
  if(x < MODEL_SIZE_X - 1 && SIGMA(lattice[x + 1][y][z]) != sigma) {
    change += bondEnergy[TYPE(lattice[x + 1][y][z])][type];
  }
  if(y > 0 && SIGMA(lattice[x][y - 1][z]) != sigma) {
    change += bondEnergy[TYPE(lattice[x][y - 1][z])][type];
  }
  if(y < MODEL_SIZE_Y - 1 && SIGMA(lattice[x][y + 1][z]) != sigma) {
    change += bondEnergy[TYPE(lattice[x][y + 1][z])][type];
  }
  if(z > 0 && SIGMA(lattice[x][y][z - 1]) != sigma) {
    change += bondEnergy[TYPE(lattice[x][y][z - 1])][type];
  }
  if(z < MODEL_SIZE_Z - 1 && SIGMA(lattice[x][y][z + 1]) != sigma) {
    change += bondEnergy[TYPE(lattice[x][y][z + 1])][type];
  }

  DEBUG(printf("\t\t\tAdhesion inser from [%i, %i, %i] sig %i: %i\n", x, y, z, sigma, change);)
  return change;
}

void calculateNextStep(int64_t ***lattice, int cellCount, cell_info_t *cells) {
  DEBUG(printf("calculating next state, global energy before: %lli", globalEnergy);)
  int i, j, x, y, z;
  site_t *neighborList = malloc(6 * sizeof(site_t));
  for(i = 0; i < cellCount; i++) {
    // reset medium
    cells[0].volume = targetVolume[0];
    cells[0].membraneArea = targetMembrane[0];

    x = rand() % MODEL_SIZE_X;
    y = rand() % MODEL_SIZE_Y;
    z = rand() % MODEL_SIZE_Z;
    //DEBUG(printf("bbb");)
    j = 0;
    if(x < MODEL_SIZE_X - 1 && SIGMA(lattice[x+1][y][z]) != SIGMA(lattice[x][y][z])) {
      neighborList[j].x = x+1; neighborList[j].y = y; neighborList[j++].z = z;
    }
    if(x > 0 && SIGMA(lattice[x-1][y][z]) != SIGMA(lattice[x][y][z])) {
      neighborList[j].x = x-1; neighborList[j].y = y; neighborList[j++].z = z;
    }
    if(y < MODEL_SIZE_Y - 1 && SIGMA(lattice[x][y+1][z]) != SIGMA(lattice[x][y][z])) {
      neighborList[j].x = x; neighborList[j].y = y+1; neighborList[j++].z = z;
    }
    if(y > 0 && SIGMA(lattice[x][y-1][z]) != SIGMA(lattice[x][y][z])) {
      neighborList[j].x = x; neighborList[j].y = y-1; neighborList[j++].z = z;
    }
    if(z < MODEL_SIZE_Z - 1 && SIGMA(lattice[x][y][z+1]) != SIGMA(lattice[x][y][z])) {
      neighborList[j].x = x; neighborList[j].y = y; neighborList[j++].z = z+1;
    }
    if(z > 0 && SIGMA(lattice[x][y][z-1]) != SIGMA(lattice[x][y][z])) {
      neighborList[j].x = x; neighborList[j].y = y; neighborList[j++].z = z-1;
    }
    //DEBUG(printf("aaa");)
    if(j > 0) {
      int choice = rand() % j;
      checkAndMoveToSite(x, y, z, neighborList[choice].x, neighborList[choice].y, neighborList[choice].z, lattice, cells);
    }
    else {
      i--;
    }
    
  }
  //free(neighborList);
  DEBUG(printf("calculating next state done, global energy after: %lli, numCells: %i", globalEnergy.total, numCells);)
}

// checks whether cell should expand to target lattie and performs the expansion if it should
// returns 1 if expansion has been performed, and 0 otherwise
int checkAndMoveToSite(int x1, int y1, int z1, int x2, int y2, int z2, int64_t ***lattice, cell_info_t *cells) {
  DEBUG(printf("\tinvading attempt start\n");)
  if(x2 >= 0 && x2 < MODEL_SIZE_X && y2 >= 0 && y2 < MODEL_SIZE_Y && z2 >= 0 && z2 < MODEL_SIZE_Z &&
      SIGMA(lattice[x2][y2][z2]) != SIGMA(lattice[x1][y1][z1])) {
    int sigmaTarget = SIGMA(lattice[x2][y2][z2]);
    int sigmaSource = SIGMA(lattice[x1][y1][z1]);
    DEBUG(printf("\tsigma %i to %i, vol %i and %i\n", sigmaSource, sigmaTarget, cells[sigmaSource].volume, cells[sigmaTarget].volume);)
    energy_t newEnergy = getNewEnergyLatticeChange(x2, y2, z2, lattice, cells, sigmaSource, sigmaTarget);

    int deltaEnergy = newEnergy.total - globalEnergy.total;
    DEBUG(printf("\tdelta energy: %i\n", deltaEnergy);)
    if(deltaEnergy < ENERGY_TRESHOLD || ((float)rand()) / ((float)RAND_MAX) < exp(-(deltaEnergy + ENERGY_TRESHOLD) / TEMPERATURE)) {
      globalEnergy.total = newEnergy.total; 
      globalEnergy.adhesion = newEnergy.adhesion;
      globalEnergy.membrane = newEnergy.membrane;
      globalEnergy.volume = newEnergy.volume;
      DEBUG(printf("\tinvading sigma %i [%i, %i, %i] to sigma %i [%i, %i, %i]\n", sigmaSource, x1, y1, z1, sigmaTarget, x2, y2, z2);)
      
      changeSiteOwner(x2, y2, z2, lattice, cells, sigmaSource, sigmaTarget);

      lattice[x2][y2][z2] = lattice[x1][y1][z1];
      DEBUG(printf("\tinvading completed, checking whether to split\n");)
      if(sigmaSource != 0 && cells[sigmaSource].volume > targetVolume[cells[sigmaSource].type] * CELL_OVERSIZE_SPLIT_TRESHOLD) {
        splitCell(x2, y2, z2, lattice, cells, sigmaSource);
      }
      return 1;
    }
  }
  DEBUG(printf("\tinvading attempt failed\n");)
  return 0;
}

void splitCell(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigma) {
  int i, subcellsToSplit = cells[sigma].volume / 2;
  DEBUG(printf("\t\tsplitting cell %i into %i, starting at [%i, %i, %i]. Subcells to split: %i\n", sigma, numCells, x, y, z, subcellsToSplit);)

  cells[numCells + 1].volume = 0;
  cells[numCells + 1].membraneArea = 0;
  cells[numCells + 1].type = cells[sigma].type;
  cells[numCells + 1].subcells = (site_t *)malloc(targetVolume[cells[numCells + 1].type] * 2 * sizeof(site_t));

  energy_t newEnergy = getNewEnergyLatticeChange(x, y, z, lattice, cells, numCells + 1, sigma);
  globalEnergy.total = newEnergy.total; 
  globalEnergy.adhesion = newEnergy.adhesion;
  globalEnergy.membrane = newEnergy.membrane;
  globalEnergy.volume = newEnergy.volume;
  numCells++;
  lattice[x][y][z] = SET_SIGMA(lattice[x][y][z], numCells);

  changeSiteOwner(x, y, z, lattice, cells, numCells, sigma);
  
  int numSubcells = 0;
  while(subcellsToSplit > 0 && numSubcells != cells[numCells].volume) {
    numSubcells = cells[numCells].volume;
    energy_t newEnergy = {0, 0, 0, 0};
    for(i = 0; i < numSubcells; i++) {
      site_t *site = &cells[numCells].subcells[i];
      if((site->x) + 1 < MODEL_SIZE_X && SIGMA(lattice[(site->x) + 1][site->y][site->z]) == sigma) {
        globalEnergy = getNewEnergyLatticeChange((site->x) + 1, site->y, site->z, lattice, cells, numCells, sigma);
        changeSiteOwner((site->x) + 1, site->y, site->z, lattice, cells, numCells, sigma);
        lattice[(site->x) + 1][site->y][site->z] = SET_SIGMA(lattice[(site->x) + 1][site->y][site->z], numCells);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->x) - 1 >= 0 && SIGMA(lattice[(site->x) - 1][site->y][site->z]) == sigma) {
        globalEnergy = getNewEnergyLatticeChange((site->x) - 1, site->y, site->z, lattice, cells, numCells, sigma);
        changeSiteOwner((site->x) - 1, site->y, site->z, lattice, cells, numCells, sigma);
        lattice[(site->x) - 1][site->y][site->z] = SET_SIGMA(lattice[(site->x) - 1][site->y][site->z], numCells);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->y) + 1 < MODEL_SIZE_Y && SIGMA(lattice[site->x][(site->y) + 1][site->z]) == sigma) {
        globalEnergy = getNewEnergyLatticeChange(site->x, (site->y) + 1, site->z, lattice, cells, numCells, sigma);
        changeSiteOwner(site->x, (site->y) + 1, site->z, lattice, cells, numCells, sigma);
        lattice[(site->x)][(site->y) + 1][site->z] = SET_SIGMA(lattice[(site->x)][(site->y) + 1][site->z], numCells);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->y) - 1 >= 0 && SIGMA(lattice[site->x][(site->y) - 1][site->z]) == sigma) {
        globalEnergy = getNewEnergyLatticeChange(site->x, (site->y) - 1, site->z, lattice, cells, numCells, sigma);
        changeSiteOwner(site->x, (site->y) - 1, site->z, lattice, cells, numCells, sigma);
        lattice[(site->x)][(site->y) - 1][site->z] = SET_SIGMA(lattice[(site->x)][(site->y) - 1][site->z], numCells);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->z) + 1 < MODEL_SIZE_Z && SIGMA(lattice[site->x][site->y][(site->z) + 1]) == sigma) {
        globalEnergy = getNewEnergyLatticeChange(site->x, site->y, (site->z) + 1, lattice, cells, numCells, sigma);
        changeSiteOwner(site->x, site->y, (site->z) + 1, lattice, cells, numCells, sigma);
        lattice[(site->x)][site->y][(site->z) + 1] = SET_SIGMA(lattice[(site->x)][site->y][(site->z) + 1], numCells);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->z) - 1 >= 0 && SIGMA(lattice[site->x][site->y][(site->z) - 1]) == sigma) {
        globalEnergy = getNewEnergyLatticeChange(site->x, site->y, (site->z) - 1, lattice, cells, numCells, sigma);
        changeSiteOwner(site->x, site->y, (site->z) - 1, lattice, cells, numCells, sigma);
        lattice[(site->x)][site->y][(site->z) - 1] = SET_SIGMA(lattice[(site->x)][site->y][(site->z) - 1], numCells);
        subcellsToSplit--;
      }
    }
  }
}

energy_t getNewEnergyLatticeChange(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigmaTo, int sigmaFrom) {
  energy_t newEnergy = {0, 0, 0, 0};
  int volumePartSource = 0, volumePartTarget = 0, membranePartTarget = 0, membranePartSource = 0;

  if(sigmaFrom != 0) {
    volumePartTarget = pow(cells[sigmaFrom].volume - targetVolume[cells[sigmaFrom].type], 2);
    membranePartTarget = pow(cells[sigmaFrom].membraneArea - targetMembrane[cells[sigmaFrom].type], 2);
  }
  if(sigmaTo != 0) {
    volumePartSource = pow(cells[sigmaTo].volume - targetVolume[cells[sigmaTo].type], 2);
    membranePartSource = pow(cells[sigmaTo].membraneArea - targetMembrane[cells[sigmaTo].type], 2);
  }
  long long volumePart = globalEnergy.volume / VOLUME_RESISTANCE - volumePartTarget - volumePartSource;
  long long membranePart =  globalEnergy.membrane / MEMBRANE_RESISTANCE - membranePartSource - membranePartTarget;

  if(sigmaTo != 0) {
    volumePart += pow(cells[sigmaTo].volume + 1 - targetVolume[cells[sigmaTo].type], 2);
    membranePart += pow(getMembraneChangeInsert(x, y, z, lattice, sigmaTo) + cells[sigmaTo].membraneArea - targetMembrane[cells[sigmaTo].type], 2);
  }
  if(sigmaFrom != 0) {
    volumePart += pow(cells[sigmaFrom].volume - 1 - targetVolume[cells[sigmaFrom].type], 2);
    membranePart += pow(getMembraneChangeRemove(x, y, z, lattice, sigmaFrom) + cells[sigmaFrom].membraneArea - targetMembrane[cells[sigmaFrom].type], 2);
  }
  volumePart = volumePart * VOLUME_RESISTANCE;
  membranePart = membranePart * MEMBRANE_RESISTANCE;

  // add energy of adhesion change
  newEnergy.adhesion = globalEnergy.adhesion;
  newEnergy.adhesion += getAdhesionChangeInsert(x, y, z, lattice, sigmaTo);
  newEnergy.adhesion += getAdhesionChangeRemove(x, y, z, lattice, sigmaFrom);

  newEnergy.membrane = membranePart;
  newEnergy.volume = volumePart;
  newEnergy.total = newEnergy.membrane + newEnergy.volume + newEnergy.adhesion;
  return newEnergy;
}

void changeSiteOwner(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigmaTo, int sigmaFrom) {
  int i;

  if(sigmaTo != 0) {
    // DEBUG(printf("\taddr %i", &cells[sigmaTo].subcells[cells[sigmaTo].volume - 1]);)
    // DEBUG(printf(", val %i\n", cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].val);)
    // move site from target cell definition to source
    cells[sigmaTo].volume++;
    cells[sigmaTo].membraneArea += getMembraneChangeInsert(x, y, z, lattice, sigmaTo);
    cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].x = x;
    cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].y = y;
    cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].z = z;
    cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].val = &lattice[x][y][z];
  }
  
  if(sigmaFrom != 0) {
    cells[sigmaFrom].volume--;
    cells[sigmaFrom].membraneArea += getMembraneChangeRemove(x, y, z, lattice, sigmaFrom);
    for(i = 0; i <= cells[sigmaFrom].volume; i++) {
      if(cells[sigmaFrom].subcells[i].val == &lattice[x][y][z]) {
        cells[sigmaFrom].subcells[i].x = cells[sigmaFrom].subcells[cells[sigmaFrom].volume].x;
        cells[sigmaFrom].subcells[i].y = cells[sigmaFrom].subcells[cells[sigmaFrom].volume].y;
        cells[sigmaFrom].subcells[i].z = cells[sigmaFrom].subcells[cells[sigmaFrom].volume].z;
        cells[sigmaFrom].subcells[i].val = cells[sigmaFrom].subcells[cells[sigmaFrom].volume].val;
      }
    }
  }
}
