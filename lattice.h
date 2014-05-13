#define SIGMA(site)               (site & 0xFFFFFFF) 
#define SET_SIGMA(c, d)           ((c & 0xFFFFFFFFF0000000) | d) 

#define TYPE(site)                ((site & 0xF0000000) >> 28)
#define MEDIUM                    0
#define VASCULAR                  1
#define TUMOR_NORM                2
#define TUMOR_NECROSIS            3
#define TUMOR_STEM                4
#define TYPE5                     5
#define SET_TYPE(c, t)            ((c & 0xFFFFFFFF0FFFFFFF) | (t << 28))   

int64_t ***initLattice();
void drawLatticeSite(int x, int y, int z, int64_t value, float xrot, float yrot, int64_t ***lattice);