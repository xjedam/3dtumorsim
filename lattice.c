#include "3dtumorsim.h"
#include "potts.h"
#include "lattice.h"

float viewScaleX = 1.0f / MODEL_SIZE_X / 2;
float viewScaleY = 1.0f / MODEL_SIZE_Y / 2;
float viewScaleZ = 1.0f / MODEL_SIZE_Z / 2;
float translationX = 1.0f / MODEL_SIZE_X;
float translationY = -1.0f / MODEL_SIZE_Y;
float translationZ = -1.0f / MODEL_SIZE_Z;
  
int64_t ***initLattice() {
  int i, j, k;
  int64_t ***ptr;
  ptr = (int64_t ***)malloc(MODEL_SIZE_X * sizeof(int64_t **));
  for(i = 0; i < MODEL_SIZE_X; i++) {
    ptr[i] = (int64_t **)malloc(MODEL_SIZE_Y * sizeof(int64_t *));
    for(j = 0; j < MODEL_SIZE_Y; j++) {
      ptr[i][j] = (int64_t *)malloc(MODEL_SIZE_Z * sizeof(int64_t));
      for(k = 0; k < MODEL_SIZE_Z; k++) {
        ptr[i][j][k] = 0;
      }
    }
  }
  
  return ptr;
}

void drawCells(float xrot, float yrot, int64_t ***lattice, cell_info_t *cells) {
  int i, j, x, y, z, sigma;
  for(i = 1; i <= numCells; i++) {
    sigma = i;
    DEBUG(printf("\t\tprinting cell sig %i with volume %i, membrane: %i, type: %i.\n", i, cells[i].volume, cells[i].membraneArea, cells[i].type);)
    for(j = 0; j < cells[i].volume; j++) {
      x = cells[i].subcells[j].x;
      y = cells[i].subcells[j].y;
      z = cells[i].subcells[j].z;
      glLoadIdentity();                       // Reset the model-view matrix
      glTranslatef((x - (MODEL_SIZE_X/2)) * translationX, (y - (MODEL_SIZE_Y/2)) * translationY, (z - (MODEL_SIZE_Z/2)) * translationZ); 
      glBegin(GL_QUADS);
        switch(cells[i].type){
          case MEDIUM:
            return;
            break;
          case VASCULAR:
            glColor4f(1.0f, 0.0f, 0.0f, 0.4f);    
            break;
          case TUMOR_NORM:
            glColor4f(0.0f, 1.0f, 0.0f, 0.4f);    
            break;
          case TUMOR_NECROSIS:
            glColor4f(0.0f, 0.0f, 1.0f, 0.4f);    
            break;
          case TUMOR_STEM:
            glColor4f(0.2f, 0.6f, 0.8f, 0.4f);    
            break;
          default:
            glColor4f(1.0f, 1.0f, 1.0f, 0.4f);
        }

        // Top face
        // DEBUG(printf("\t\t\ty+1 is sig: %lli, type: %lli", SIGMA(lattice[x][y + 1][z]), TYPE(lattice[x][y + 1][z]));)
        if(y >= MODEL_SIZE_Y - 1 || SIGMA(lattice[x][y + 1][z]) != sigma) {
          glVertex3f( viewScaleX, viewScaleY, -viewScaleZ);
          glVertex3f(-viewScaleX, viewScaleY, -viewScaleZ);
          glVertex3f(-viewScaleX, viewScaleY,  viewScaleZ);
          glVertex3f( viewScaleX, viewScaleY,  viewScaleZ);
        }

        // Bottom face 
        // DEBUG(printf("\t\t\ty-1 is sig: %lli, type: %lli", SIGMA(lattice[x][y - 1][z]), TYPE(lattice[x][y + 1][z]));)
        if(y == 0 || SIGMA(lattice[x][y - 1][z]) != sigma) { 
          glVertex3f( viewScaleX, -viewScaleY,  viewScaleZ);
          glVertex3f(-viewScaleX, -viewScaleY,  viewScaleZ);
          glVertex3f(-viewScaleX, -viewScaleY, -viewScaleZ);
          glVertex3f( viewScaleX, -viewScaleY, -viewScaleZ);
        }

        // Front face 
        // DEBUG(printf("\t\t\tz+1 is sig: %lli, type: %lli", SIGMA(lattice[x][y][z + 1]), TYPE(lattice[x][y][z + 1]));)
        if(z >= MODEL_SIZE_Z - 1 || SIGMA(lattice[x][y][z + 1]) != sigma) {
          glVertex3f( viewScaleX,  viewScaleY, viewScaleZ);
          glVertex3f(-viewScaleX,  viewScaleY, viewScaleZ);
          glVertex3f(-viewScaleX, -viewScaleY, viewScaleZ);
          glVertex3f( viewScaleX, -viewScaleY, viewScaleZ);
        }

        // Back face 
        // DEBUG(printf("\t\t\tz-1 is sig: %lli, type: %lli", SIGMA(lattice[x][y][z - 1]), TYPE(lattice[x][y][z - 1]));)
        if(z == 0 || SIGMA(lattice[x][y][z - 1]) != sigma) {
          glVertex3f( viewScaleX, -viewScaleY, -viewScaleZ);
          glVertex3f(-viewScaleX, -viewScaleY, -viewScaleZ);
          glVertex3f(-viewScaleX,  viewScaleY, -viewScaleZ);
          glVertex3f( viewScaleX,  viewScaleY, -viewScaleZ);
        }

        // Left face 
        // DEBUG(printf("\t\t\tx-1 is sig: %lli, type: %lli", SIGMA(lattice[x - 1][y][z]), TYPE(lattice[x - 1][y][z]));)
        if(x == 0 || SIGMA(lattice[x - 1][y][z]) != sigma) {
          glVertex3f(-viewScaleX,  viewScaleY,  viewScaleZ);
          glVertex3f(-viewScaleX,  viewScaleY, -viewScaleZ);
          glVertex3f(-viewScaleX, -viewScaleY, -viewScaleZ);
          glVertex3f(-viewScaleX, -viewScaleY,  viewScaleZ);
        }

        // Right face 
        // DEBUG(printf("\t\t\tx+1 is sig: %lli, type: %lli", SIGMA(lattice[x + 1][y][z]), TYPE(lattice[x + 1][y][z]));)
        if(x >= MODEL_SIZE_X - 1 || SIGMA(lattice[x + 1][y][z]) != sigma) {
          glVertex3f(viewScaleX,  viewScaleY, -viewScaleZ);
          glVertex3f(viewScaleX,  viewScaleY,  viewScaleZ);
          glVertex3f(viewScaleX, -viewScaleY,  viewScaleZ);
          glVertex3f(viewScaleX, -viewScaleY, -viewScaleZ);
        }
        
      glEnd();
    }
  }
}

void drawLatticeSite(int x, int y, int z, int64_t value, float xrot, float yrot, int64_t ***lattice) {
  glLoadIdentity();                       // Reset the model-view matrix
  glTranslatef((x - (MODEL_SIZE_X/2)) * translationX, (y - (MODEL_SIZE_Y/2)) * translationY, (z - (MODEL_SIZE_Z/2)) * translationZ); 
  
  glBegin(GL_QUADS);
    switch(TYPE(value)){
      case MEDIUM:
        return;
        break;
      case VASCULAR:
        glColor4f(1.0f, 0.0f, 0.0f, 1.0f);    
        break;
      case TUMOR_NORM:
        glColor4f(0.0f, 1.0f, 0.0f, 1.0f);    
        break;
      case TUMOR_NECROSIS:
        glColor4f(0.0f, 0.0f, 1.0f, 1.0f);    
        break;
      case TUMOR_STEM:
        glColor4f(0.2f, 0.6f, 0.8f, 1.0f);    
        break;
      default:
        glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    }
    int sigma = SIGMA(lattice[x][y][z]);

    // Top face
    if(y == MODEL_SIZE_Y - 1 || SIGMA(lattice[x][y + 1][z]) != sigma) {
      glVertex3f( viewScaleX, viewScaleY, -viewScaleZ);
      glVertex3f(-viewScaleX, viewScaleY, -viewScaleZ);
      glVertex3f(-viewScaleX, viewScaleY,  viewScaleZ);
      glVertex3f( viewScaleX, viewScaleY,  viewScaleZ);
    }

    // Bottom face 
    if(y == 0 || SIGMA(lattice[x][y - 1][z]) != sigma) { 
      glVertex3f( viewScaleX, -viewScaleY,  viewScaleZ);
      glVertex3f(-viewScaleX, -viewScaleY,  viewScaleZ);
      glVertex3f(-viewScaleX, -viewScaleY, -viewScaleZ);
      glVertex3f( viewScaleX, -viewScaleY, -viewScaleZ);
    }

    // Front face 
    if(z == MODEL_SIZE_Z - 1 || SIGMA(lattice[x][y][z + 1]) != sigma) {
      glVertex3f( viewScaleX,  viewScaleY, viewScaleZ);
      glVertex3f(-viewScaleX,  viewScaleY, viewScaleZ);
      glVertex3f(-viewScaleX, -viewScaleY, viewScaleZ);
      glVertex3f( viewScaleX, -viewScaleY, viewScaleZ);
    }

    // Back face 
    if(z == 0 || SIGMA(lattice[x][y][z - 1]) != sigma) {
      glVertex3f( viewScaleX, -viewScaleY, -viewScaleZ);
      glVertex3f(-viewScaleX, -viewScaleY, -viewScaleZ);
      glVertex3f(-viewScaleX,  viewScaleY, -viewScaleZ);
      glVertex3f( viewScaleX,  viewScaleY, -viewScaleZ);
    }

    // Left face 
    if(x == 0 || SIGMA(lattice[x - 1][y][z]) != sigma) {
      glVertex3f(-viewScaleX,  viewScaleY,  viewScaleZ);
      glVertex3f(-viewScaleX,  viewScaleY, -viewScaleZ);
      glVertex3f(-viewScaleX, -viewScaleY, -viewScaleZ);
      glVertex3f(-viewScaleX, -viewScaleY,  viewScaleZ);
    }

    // Right face 
    if(x == MODEL_SIZE_X - 1 || SIGMA(lattice[x + 1][y][z]) != sigma) {
      glVertex3f(viewScaleX,  viewScaleY, -viewScaleZ);
      glVertex3f(viewScaleX,  viewScaleY,  viewScaleZ);
      glVertex3f(viewScaleX, -viewScaleY,  viewScaleZ);
      glVertex3f(viewScaleX, -viewScaleY, -viewScaleZ);
    }
    
  glEnd();
}