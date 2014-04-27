#include "3dtumorsim.h"
#include "lattice.h"
#include "potts.h"

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
    cells[sigma].membraneArea += getMembraneChange(x, y, z, lattice);
    if(sigma > max) {
      max = sigma;
    }
	}

	fclose(fp);
  return max;
}

// gets the change of the membrane area when x, y, z was inserted
int getMembraneChange(int x, int y, int z, int64_t ***lattice) {
  int change = -6;
  int sigma = SIGMA(lattice[x][y][z]);

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