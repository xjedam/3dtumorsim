#include "3dtumorsim.h"
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