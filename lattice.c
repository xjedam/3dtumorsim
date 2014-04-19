#include "3dtumorsim.h"
#include "lattice.h"

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

void drawLatticeSite(int x, int y, int z, int64_t value) {
  
  glLoadIdentity();                       // Reset the model-view matrix
  glTranslatef(x * SITE_SIZE, y * SITE_SIZE, z * SITE_SIZE);                  // Move left and into the screen
  
  glBegin(GL_QUADS);
    switch(TYPE(value)){
      case MEDIUM:
        return;
        break;
      case VASCULAR:
        glColor4f(1.0f, 0.0f, 0.0f, 0.4f);    
        break;
      case TUMOR_NORM:
        glColor4f(0.0f, 1.0f, 0.0f, 0.4f);    
        break;
      default:
        glColor4f(1.0f, 1.0f, 1.0f, 0.4f);
    }

    glVertex3f( 1.0f, 1.0f, -1.0f);
    glVertex3f(-1.0f, 1.0f, -1.0f);
    glVertex3f(-1.0f, 1.0f,  1.0f);
    glVertex3f( 1.0f, 1.0f,  1.0f);

    // Bottom face 
    glVertex3f( 1.0f, -1.0f,  1.0f);
    glVertex3f(-1.0f, -1.0f,  1.0f);
    glVertex3f(-1.0f, -1.0f, -1.0f);
    glVertex3f( 1.0f, -1.0f, -1.0f);

    // Front face  
    glVertex3f( 1.0f,  1.0f, 1.0f);
    glVertex3f(-1.0f,  1.0f, 1.0f);
    glVertex3f(-1.0f, -1.0f, 1.0f);
    glVertex3f( 1.0f, -1.0f, 1.0f);

    // Back face 
    glVertex3f( 1.0f, -1.0f, -1.0f);
    glVertex3f(-1.0f, -1.0f, -1.0f);
    glVertex3f(-1.0f,  1.0f, -1.0f);
    glVertex3f( 1.0f,  1.0f, -1.0f);

    // Left face 
    glVertex3f(-1.0f,  1.0f,  1.0f);
    glVertex3f(-1.0f,  1.0f, -1.0f);
    glVertex3f(-1.0f, -1.0f, -1.0f);
    glVertex3f(-1.0f, -1.0f,  1.0f);

    // Right face 
    glVertex3f(1.0f,  1.0f, -1.0f);
    glVertex3f(1.0f,  1.0f,  1.0f);
    glVertex3f(1.0f, -1.0f,  1.0f);
    glVertex3f(1.0f, -1.0f, -1.0f);
    
  glEnd();
}