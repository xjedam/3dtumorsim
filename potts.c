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
		cells[i].type = 0;
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
		lattice[x][y][z] = SET_TYPE(SET_SIGMA(0,sigma), type);
		cells[sigma].temperature = temperature;
		cells[sigma].type = type;
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
  long long newEnergy = globalEnergy;

  if(x2 >= 0 && x2 < MODEL_SIZE_X && y2 >= 0 && y2 < MODEL_SIZE_Y && z2 >= 0 && z2 < MODEL_SIZE_Z &&
      SIGMA(lattice[x2][y2][z2]) != SIGMA(lattice[x1][y1][z1])) {
    int sigmaTarget = SIGMA(lattice[x2][y2][z2]);
    int sigmaSource = SIGMA(lattice[x1][y1][z1]);

    int membranePartTarget = pow(cells[sigmaTarget].membraneArea - targetMembrane[cells[sigmaTarget].type], 2);
    int membranePartSource = pow(cells[sigmaSource].membraneArea - targetMembrane[cells[sigmaSource].type], 2);
    int membranePart =  membranePartSource + membranePartTarget;
    int volumePartTarget = pow(cells[sigmaTarget].volume - targetVolume[cells[sigmaTarget].type], 2);
    int volumePartSource = pow(cells[sigmaSource].volume - targetVolume[cells[sigmaSource].type], 2);
    int volumePart = volumePartTarget + volumePartSource;

    int oldEnergyValues = MEMBRANE_RESISTANCE * membranePart + VOLUME_RESISTANCE * volumePart;

    // subtracts old cell membrane and volume values
    newEnergy -= oldEnergyValues;
    // add energy of new source cell and target cell membrane state
    newEnergy += (getMembraneChangeInsert(x2, y2, z2, lattice, sigmaSource) + membranePartSource) * MEMBRANE_RESISTANCE;
    newEnergy += (getMembraneChangeRemove(x2, y2, z2, lattice, sigmaTarget) + membranePartTarget) * MEMBRANE_RESISTANCE;
    // add energy of new source cell and target cell volume state
    newEnergy += (volumePartSource + 1) * VOLUME_RESISTANCE;
    newEnergy += (volumePartTarget - 1) * VOLUME_RESISTANCE;
    // add energy of adhesion change
    newEnergy += getAdhesionChangeInsert(x2, y2, z2, lattice, sigmaSource);
    newEnergy += getAdhesionChangeRemove(x2, y2, z2, lattice, sigmaTarget);

    int deltaEnergy = globalEnergy - newEnergy;
    if(deltaEnergy < ENERGY_TRESHOLD || ((float)rand()) / ((float)RAND_MAX) < exp(-(deltaEnergy + ENERGY_TRESHOLD) / TEMPERATURE)) {
      globalEnergy = newEnergy;
      cells[sigmaSource].volume++;
      cells[sigmaTarget].volume--;
      cells[sigmaSource].membraneArea += getMembraneChangeInsert(x2, y2, z2, lattice, sigmaSource);
      cells[sigmaTarget].membraneArea -= getMembraneChangeRemove(x2, y2, z2, lattice, sigmaTarget);
      lattice[x2][y2][z2] = lattice[x1][y1][z1];
      return 1;
    }

  }
  return 0;
}

void splitCell(int x, int y, int z, int64_t ***lattice, cell_info_t *cells) {
  
}