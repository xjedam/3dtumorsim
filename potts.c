#include "3dtumorsim.h"
#include "lattice.h"
#include "potts.h"
#include <math.h>

// defining bond energy
int bondEnergy[6][6] = {
  {0, 12, 10, 15, 10, 0},       // MEDIUM
  {12, 5, 30, 30, 30, 0},       // VASCULAR
  {10, 30, 8, 8, 8, 0},         // TUMOR_NORM
  {15, 30, 8, 3, 8, 0},         // TUMOR_NECROSIS
  {10, 30, 8, 8, 8, 0},         // TUMOR_STEM
  {0, 0, 0, 0, 0, 0},
};

// target membrane area definition
int targetMembrane[6] = {
  30,                           // MEDIUM
  40,                           // VASCULAR
  30,                           // TUMOR_NORM
  30,                           // TUMOR_NECROSIS
  30,                           // TUMOR_STEM
  0
};

// target volume definition
int targetVolume[6] = {
  15,                           // MEDIUM
  24,                           // VASCULAR
  15,                           // TUMOR_NORM
  15,                           // TUMOR_NECROSIS
  15,                           // TUMOR_STEM
  0
};

long long globalEnergy = 0;

cell_info_t *initCells() {
	int i, size = MODEL_SIZE_X * MODEL_SIZE_Y * MODEL_SIZE_Z / 8;
	cell_info_t *cells = (cell_info_t *)malloc(size * sizeof(cell_info_t));
	for(i = 0; i < size; i++) {
		cells[i].temperature = 36.7;
		cells[i].type = MEDIUM;
		cells[i].volume = 30;
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
	int x, y, z, sigma, temperature, type, i, max = 0, j, k;
	while(fscanf(fp, "%i%i%i:%i%i%f", &x, &y, &z, &type, &sigma, &temperature) != EOF) {
    if(cells[sigma].type != type) {
      cells[sigma].subcells = (site_t *)malloc(targetVolume[type] * 2 * sizeof(site_t));
      cells[sigma].type = type;
      cells[sigma].temperature = temperature;
      cells[sigma].volume = 0;
      cells[sigma].membraneArea = 0;
    }
		lattice[x][y][z] = SET_TYPE(SET_SIGMA(0,sigma), type);
    cells[sigma].subcells[cells[sigma].volume].val = &lattice[x][y][z];
    cells[sigma].subcells[cells[sigma].volume].x = x;
    cells[sigma].subcells[cells[sigma].volume].y = y;
    cells[sigma].subcells[cells[sigma].volume].z = z;
		cells[sigma].volume += 1;
    cells[sigma].membraneArea += getMembraneChangeInsert(x, y, z, lattice, sigma);
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
  globalEnergy = MEMBRANE_RESISTANCE * membranePart + VOLUME_RESISTANCE * volumePart;

  // add global energy: adhesion part
  for(i = 0; i < MODEL_SIZE_X; i++) {
    for(j = 0; j < MODEL_SIZE_Y; j++) {
      for(k = 0; k < MODEL_SIZE_Z; k++) {
        if(i + 1 < MODEL_SIZE_X && SIGMA(lattice[i][j][k]) != SIGMA(lattice[i + 1][j][k])) {
          globalEnergy += bondEnergy[SIGMA(lattice[i][j][k])][SIGMA(lattice[i + 1][j][k])];
        }
        if(j + 1 < MODEL_SIZE_Y && SIGMA(lattice[i][j][k]) != SIGMA(lattice[i][j + 1][k])) {
          globalEnergy += bondEnergy[SIGMA(lattice[i][j][k])][SIGMA(lattice[i][j + 1][k])];
        }
        if(k + 1 < MODEL_SIZE_Z && SIGMA(lattice[i][j][k]) != SIGMA(lattice[i][j][k + 1])) {
          globalEnergy += bondEnergy[SIGMA(lattice[i][j][k])][SIGMA(lattice[i][j][k + 1])];
        }
      }
    }
  }

	fclose(fp);
  DEBUG(puts("loaded model");)
  return max;
}

// gets the change of the membrane area when x, y, z was inserted
int getMembraneChangeInsert(int x, int y, int z, int64_t ***lattice, int sigma) {
  int change = -6;

  if(x > 0 && SIGMA(lattice[x - 1][y][z]) != sigma) {
    change += 2;
  }
  if(x < MODEL_SIZE_X - 1 && SIGMA(lattice[x + 1][y][z]) != sigma) {
    change += 2;
  }
  if(y > 0 && SIGMA(lattice[x][y - 1][z]) != sigma) {
    change += 2;
  }
  if(y < MODEL_SIZE_Y - 1 && SIGMA(lattice[x][y + 1][z]) != sigma) {
    change += 2;
  }
  if(z > 0 && SIGMA(lattice[x][y][z - 1]) != sigma) {
    change += 2;
  }
  if(z < MODEL_SIZE_Z - 1 && SIGMA(lattice[x][y][z + 1]) != sigma) {
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
  int change = 0;

  if(x > 0 && SIGMA(lattice[x - 1][y][z]) != sigma) {
    change -= bondEnergy[SIGMA(lattice[x - 1][y][z])][sigma];
  }
  if(x < MODEL_SIZE_X - 1 && SIGMA(lattice[x + 1][y][z]) != sigma) {
    change -= bondEnergy[SIGMA(lattice[x + 1][y][z])][sigma];
  }
  if(y > 0 && SIGMA(lattice[x][y - 1][z]) != sigma) {
    change -= bondEnergy[SIGMA(lattice[x][y - 1][z])][sigma];
  }
  if(y < MODEL_SIZE_Y - 1 && SIGMA(lattice[x][y + 1][z]) != sigma) {
    change -= bondEnergy[SIGMA(lattice[x][y + 1][z])][sigma];
  }
  if(z > 0 && SIGMA(lattice[x][y][z - 1]) != sigma) {
    change -= bondEnergy[SIGMA(lattice[x][y][z - 1])][sigma];
  }
  if(z < MODEL_SIZE_Z - 1 && SIGMA(lattice[x][y][z + 1]) != sigma) {
    change -= bondEnergy[SIGMA(lattice[x][y][z + 1])][sigma];
  }

  return change;
}

int getAdhesionChangeInsert(int x, int y, int z, int64_t ***lattice, int sigma) {
  int change = 0;

  if(x > 0 && SIGMA(lattice[x - 1][y][z]) != sigma) {
    change += bondEnergy[SIGMA(lattice[x - 1][y][z])][sigma];
  }
  if(x < MODEL_SIZE_X - 1 && SIGMA(lattice[x + 1][y][z]) != sigma) {
    change += bondEnergy[SIGMA(lattice[x + 1][y][z])][sigma];
  }
  if(y > 0 && SIGMA(lattice[x][y - 1][z]) != sigma) {
    change += bondEnergy[SIGMA(lattice[x][y - 1][z])][sigma];
  }
  if(y < MODEL_SIZE_Y - 1 && SIGMA(lattice[x][y + 1][z]) != sigma) {
    change += bondEnergy[SIGMA(lattice[x][y + 1][z])][sigma];
  }
  if(z > 0 && SIGMA(lattice[x][y][z - 1]) != sigma) {
    change += bondEnergy[SIGMA(lattice[x][y][z - 1])][sigma];
  }
  if(z < MODEL_SIZE_Z - 1 && SIGMA(lattice[x][y][z + 1]) != sigma) {
    change += bondEnergy[SIGMA(lattice[x][y][z + 1])][sigma];
  }

  return change;
}

void calculateNextStep(int64_t ***lattice, int cellCount, cell_info_t *cells) {
  DEBUG(puts("calculating next state");)
  int i, j, x, y, z, xModiff, yModiff, zModiff;
  for(i = 0; i < cellCount; i++) {
    // reset medium
    cells[0].volume = targetVolume[0];
    cells[0].membraneArea = targetMembrane[0];

    x = rand() % MODEL_SIZE_X;
    y = rand() % MODEL_SIZE_Y;
    z = rand() % MODEL_SIZE_Z;

    xModiff = (rand() % 2) * 2 - 1;
    yModiff = (rand() % 2) * 2 - 1;
    zModiff = (rand() % 2) * 2 - 1;

    if(checkAndMoveToSite(x, y, z, x + xModiff, y, z, lattice, cells)){ continue; }
    if(checkAndMoveToSite(x, y, z, x - xModiff, y, z, lattice, cells)){ continue; }
    if(checkAndMoveToSite(x, y, z, x, y + yModiff, z, lattice, cells)){ continue; }
    if(checkAndMoveToSite(x, y, z, x, y - yModiff, z, lattice, cells)){ continue; }
    if(checkAndMoveToSite(x, y, z, x, y, z + zModiff, lattice, cells)){ continue; }
    if(checkAndMoveToSite(x, y, z, x, y, z - zModiff, lattice, cells)){ continue; }

    i--;
  }
}

// checks whether cell should expand to target lattie and performs the expansion if it should
// returns 1 if expansion has been performed, and 0 otherwise
int checkAndMoveToSite(int x1, int y1, int z1, int x2, int y2, int z2, int64_t ***lattice, cell_info_t *cells) {
  int i;  

  if(x2 >= 0 && x2 < MODEL_SIZE_X && y2 >= 0 && y2 < MODEL_SIZE_Y && z2 >= 0 && z2 < MODEL_SIZE_Z &&
      SIGMA(lattice[x2][y2][z2]) != SIGMA(lattice[x1][y1][z1])) {
    int sigmaTarget = SIGMA(lattice[x2][y2][z2]);
    int sigmaSource = SIGMA(lattice[x1][y1][z1]);

    long long newEnergy = getNewEnergyLatticeChange(x2, y2, z2, lattice, cells, sigmaSource, sigmaTarget);

    int deltaEnergy = globalEnergy - newEnergy;
    if(deltaEnergy < ENERGY_TRESHOLD || ((float)rand()) / ((float)RAND_MAX) < exp(-(deltaEnergy + ENERGY_TRESHOLD) / TEMPERATURE)) {
      globalEnergy = newEnergy;
      DEBUG(printf("\tinvading sigma %i [%i, %i, %i] to sigma %i [%i, %i, %i]\n", sigmaSource, x1, y1, z1, sigmaTarget, x2, y2, z2);)
      
      changeSiteOwner(x2, y2, z2, lattice, cells, sigmaSource, sigmaTarget);

      lattice[x2][y2][z2] = lattice[x1][y1][z1];
      DEBUG(printf("\tinvading completed\n");)
      return 1;
    }

  }
  return 0;
}

void splitCell(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigma) {
  int i, subcellsToSplit = cells[sigma].volume / 2;
  DEBUG(printf("\tsplitting cell %i into %i, starting at [%i, %i, %i]. Subcells to split: %i\n", sigma, numCells, x, y, z, subcellsToSplit);)

  cells[numCells + 1].volume = 0;
  cells[numCells + 1].membraneArea = 0;
  cells[numCells + 1].type = cells[sigma].type;
  cells[numCells + 1].subcells = (site_t *)malloc(targetVolume[cells[numCells + 1].type] * 2 * sizeof(site_t));

  long long newEnergy = getNewEnergyLatticeChange(x, y, z, lattice, cells, numCells + 1, sigma);
  globalEnergy = newEnergy;
  lattice[x][y][z] = SET_SIGMA(lattice[x][y][z], ++numCells);

  changeSiteOwner(x, y, z, lattice, cells, numCells, sigma);
  
  int numSubcells = 0;
  while(subcellsToSplit > 0 && numSubcells != cells[numCells].volume) {
    numSubcells = cells[numCells].volume;
    for(i = 0; i < numSubcells; i++) {
      site_t *site = &cells[numCells].subcells[i];
      if((site->x) + 1 < MODEL_SIZE_X && SIGMA(lattice[(site->x) + 1][site->y][site->z]) == sigma) {
        changeSiteOwner((site->x) + 1, site->y, site->z, lattice, cells, numCells, sigma);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->x) - 1 >= 0 && SIGMA(lattice[(site->x) - 1][site->y][site->z]) == sigma) {
        changeSiteOwner((site->x) - 1, site->y, site->z, lattice, cells, numCells, sigma);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->y) + 1 < MODEL_SIZE_Y && SIGMA(lattice[site->x][(site->y) + 1][site->z]) == sigma) {
        changeSiteOwner(site->x, (site->y) + 1, site->z, lattice, cells, numCells, sigma);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->y) - 1 >= 0 && SIGMA(lattice[site->x][(site->y) - 1][site->z]) == sigma) {
        changeSiteOwner(site->x, (site->y) - 1, site->z, lattice, cells, numCells, sigma);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->z) + 1 < MODEL_SIZE_Z && SIGMA(lattice[site->x][site->y][(site->z) + 1]) == sigma) {
        changeSiteOwner(site->x, site->y, (site->z) + 1, lattice, cells, numCells, sigma);
        subcellsToSplit--;
      }
      if(subcellsToSplit > 0 && (site->z) - 1 >= 0 && SIGMA(lattice[site->x][site->y][(site->z) - 1]) == sigma) {
        changeSiteOwner(site->x, site->y, (site->z) - 1, lattice, cells, numCells, sigma);
        subcellsToSplit--;
      }
    }
  }
}

long long getNewEnergyLatticeChange(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigmaTo, int sigmaFrom) {
  long long newEnergy = globalEnergy;

  int membranePartTarget = pow(cells[sigmaFrom].membraneArea - targetMembrane[cells[sigmaFrom].type], 2);
  int membranePartSource = pow(cells[sigmaTo].membraneArea - targetMembrane[cells[sigmaTo].type], 2);
  int membranePart =  membranePartSource + membranePartTarget;
  int volumePartTarget = pow(cells[sigmaFrom].volume - targetVolume[cells[sigmaFrom].type], 2);
  int volumePartSource = pow(cells[sigmaTo].volume - targetVolume[cells[sigmaTo].type], 2);
  int volumePart = volumePartTarget + volumePartSource;

  int oldEnergyValues = MEMBRANE_RESISTANCE * membranePart + VOLUME_RESISTANCE * volumePart;

  // subtracts old cell membrane and volume values
  newEnergy -= oldEnergyValues;
  // add energy of new source cell and target cell membrane state
  newEnergy += (getMembraneChangeInsert(x, y, z, lattice, sigmaTo) + membranePartSource) * MEMBRANE_RESISTANCE;
  newEnergy += (getMembraneChangeRemove(x, y, z, lattice, sigmaFrom) + membranePartTarget) * MEMBRANE_RESISTANCE;
  // add energy of new source cell and target cell volume state
  newEnergy += (volumePartSource + 1) * VOLUME_RESISTANCE;
  newEnergy += (volumePartTarget - 1) * VOLUME_RESISTANCE;
  // add energy of adhesion change
  newEnergy += getAdhesionChangeInsert(x, y, z, lattice, sigmaTo);
  newEnergy += getAdhesionChangeRemove(x, y, z, lattice, sigmaFrom);

  return newEnergy;
}

void changeSiteOwner(int x, int y, int z, int64_t ***lattice, cell_info_t *cells, int sigmaTo, int sigmaFrom) {
  int i;
  cells[sigmaTo].volume++;
  cells[sigmaFrom].volume--;
  cells[sigmaTo].membraneArea += getMembraneChangeInsert(x, y, z, lattice, sigmaTo);
  cells[sigmaFrom].membraneArea -= getMembraneChangeRemove(x, y, z, lattice, sigmaFrom);

  if(sigmaTo != 0) {
    // DEBUG(printf("\taddr %i", &cells[sigmaTo].subcells[cells[sigmaTo].volume - 1]);)
    // DEBUG(printf(", val %i\n", cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].val);)
    // move site from target cell definition to source
    cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].x = x;
    cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].y = y;
    cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].z = z;
    cells[sigmaTo].subcells[cells[sigmaTo].volume - 1].val = &lattice[x][y][z];
  }
  
  if(sigmaFrom != 0) {
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