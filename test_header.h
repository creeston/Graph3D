#include <SDL/SDL.h>

#define max2(x,y) ((x) > (y) ? (x) : (y))
#define min2(x,y) ((x) < (y) ? (x) : (y))
#define max3(x,y,z) ((x) > (y) ? max2(x,z) : max2(y,z))
#define min3(x,y,z) ((x) < (y) ? min2(x,z) : min2(y,z))
#define xwhole(x) ((int)(((x) - Xmin) / deltaX))
#define ywhole(y) ((int)(((y) - Ymin) / deltaY))
#define xreal(i) (Xmin + ((i)) * deltaX)

#define M -1000000.0
#define NVERTEX 3000
#define NTRIANGLE 2000
#define NSCREEN 40
#define BIG 1e30
#define NPOLY 400
#define NNTRSET 200


struct ver {
    double x,y,z;
    int *connect;
};

struct tr{
    int A,B,C;
    double a,b,c,h;
};

struct node {
    int jtr;
    struct node *next;
};

struct scr{
    int tr_cov;
    double tr_dist;
    struct node *start;
};

/* SDL.c*/
struct color{
    int r,g,b;
};

/////////////////////////////////////////////////////////
void linesegment(double xP, double yP, double zP, double xQ, double yQ, double zQ, int k0);
void add_linesegment(int P, int Q);
int counter_clock(int i0, int i1, int i2, double *pdist, int code);
int reflo(double *px);
int reint(int *pi);
void skipbl();
void error(char *string);
void coeff(double rho, double theta, double phi);
void viewing(double x, double y, double z, double *pxe, double *pye, double *pze);
void init_viewport(double *pXMIN, double *pXMAX, double *pYMIN, double *pYMAX0);
void opengraph();
void check();
void drawtr();

/* SDL.c */
void putpix(int x, int y, Uint32 pixel);
void init(int width, int height, int depth);
void drawpix(int x, int y, struct color C);
void clearscreen();
