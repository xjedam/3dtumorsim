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
  0,                            // MEDIUM
  40,                           // VASCULAR
  30,                           // TUMOR_NORM
  30,                           // TUMOR_NECROSIS
  30,                           // TUMOR_STEM
  0
};

// target volume definition
int targetVolume[6] = {
  0,                            // MEDIUM
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
		cells[i].volume = 0;
		cells[i].membraneArea = 0;
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
	int x, y, z, sigma, temperature, type, i, max = 0;
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

  // calculate global energy
  int membranePart = 0;
  int volumePart = 0;
  for(i = 1; i <= max; i++) {
    membranePart += pow(cells[i].membraneArea - targetMembrane[cells[i].type], 2);
    volumePart += pow(cells[i].volume - targetVolume[cells[i].type], 2);
  }
  globalEnergy = MEMBRANE_RESISTANCE * membranePart + VOLUME_RESISTANCE * volumePart;

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
    change++;
  }
  if(x < MODEL_SIZE_X - 1 && SIGMA(lattice[x + 1][y][z]) == sigma) {
    change++;
  }
  if(y > 0 && SIGMA(lattice[x][y - 1][z]) == sigma) {
    change++;
  }
  if(y < MODEL_SIZE_Y - 1 && SIGMA(lattice[x][y + 1][z]) == sigma) {
    change++;
  }
  if(z > 0 && SIGMA(lattice[x][y][z - 1]) == sigma) {
    change++;
  }
  if(z < MODEL_SIZE_Z - 1 && SIGMA(lattice[x][y][z + 1]) == sigma) {
    change++;
  }

  return change;
}

void calculateNextStep(int64_t ***lattice, int cellCount) {
  int i, x, y, z, xStart, yStart, zStart;
  for(i = 0; i < cellCount; i++) {
    x = rand() % MODEL_SIZE_X;
    y = rand() % MODEL_SIZE_Y;
    z = rand() % MODEL_SIZE_Z;
    if(TYPE(lattice[x][y][z]) != MEDIUM) {
      xStart = (rand() % 2) * 2 - 1;
      yStart = (rand() % 2) * 2 - 1;
      zStart = (rand() % 2) * 2 - 1;

      if(x + xStart >= 0 && x + xStart < MODEL_SIZE_X &&
          y + yStart >= 0 && y + yStart < MODEL_SIZE_Y &&
          z + zStart >= 0 && z + zStart < MODEL_SIZE_Z) {
        // TODO
      }
    } else {
      i--;
    }
  }
}