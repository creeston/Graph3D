#include "test_header.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <X11/Xlib.h>


extern int ntr, iaux, ipixmin, ipixmax, ipixleft, ipixright, ipix,
    jpix, jtop, jbot, j_old, l, jl, topcode[3], poly[NPOLY],
    npoly, lower[NSCREEN], upper[NSCREEN],
    low[NSCREEN], up[NSCREEN], trset[NNTRSET], ntrset;
extern double v11, v12, v13, v21, v22, v23, v32, v33, v43, d, c1, c2,
       Xrange, Yrange, Xvp_range, Yvp_range, Xmin, Xmax, Ymin, Ymax,
       deltaX, deltaY, denom, slope, Xleft[3], Xright[3], Yleft[3], Yright[3];

extern struct ver vertex[NVERTEX], *pvertex;
extern struct tr triangle[NTRIANGLE], *ptriangle;
extern struct node *pnode;
extern struct scr ScreeN[NSCREEN][NSCREEN], *pointer;



void add_linesegment(int P, int Q)
{
    int iaux, *ptr, ii, n;
    if (P > Q) {
        iaux = P;
        P = Q;
        Q = iaux;
    }

    ptr = vertex[P].connect; n = *ptr;
    for (ii = 1; ii <= n; ++ii)
        if (*(ptr + ii) == Q) return;
    ++n;
    vertex[P].connect = ptr = (int *)realloc(ptr, (n + 1) * sizeof(int));
    if (ptr == NULL)
        error("Memory alloc error \n");
    *(ptr + n) = Q; *ptr = n;
}


int counter_clock(int i0, int i1, int i2, double *pdist, int code)
{
    /*
        code = 0 : compute orientation
        code = 1 : compute a,b,c,h, store the first triangle
        code > 1 : check if next triangle is coplanar; store it
    */

    int A = abs(poly[i0]), B = abs(poly[i1]), C = abs(poly[i2]);
    double   xA, yA, zA, xB, yB, zB, xC, yC, zC, r, xdist, ydist, zdist,
             XA, YA, XB, YB, XC, YC, h0,
             DA, DB, DC, D, DAB, DAC, DBC, aux, dist, xR, yR;
    static double a,b,c,h;

    pvertex = vertex + A;
    xA = pvertex->x; yA = pvertex->y; zA = pvertex->z;
    pvertex = vertex + B;
    xB = pvertex->x; yB = pvertex->y; zB = pvertex->z;
    pvertex = vertex + C;
    xC = pvertex->x; yC = pvertex->y; zC = pvertex->z;

    h0 = xA * (yB * zC - yC * zB) - xB * (yA * zC - yC * zA) + xC * (yA * zB - yB * zA);

    if (code == 0)
        if (h0 > eps) {
            xdist = xC - xA; ydist = yC - yA; zdist = zC - zA;
            *pdist = sqrt(xdist * xdist + ydist * ydist + zdist * zdist); //!!!!!!!! -sqrt();
            return 1;
        } else return 0;

    if (code == 1) {
	a = yA * (zB - zC) - yB * (zA - zC) + yC * (zA - zB);
        b = -(xA * (zB - zC) - xB * (zA - zC) + xC * (zA - zB));
        c =  xA * (yB - yC) - xB * (yA - yC) + xC * (yA - yB);
        r = sqrt(a*a + b*b + c*c); if (r == 0) r = eps;//r = r == 0 ? eps : r;
        a /= r; b /= r; c /= r; h = h0 / r;
    } else if (fabs(a * xC + b * yC + c * zC - h) > 0.001 * fabs(h))
        error("Incorrectly specified polygon \n");

    if (ntr == NTRIANGLE) error("Too many triangles");

    ptriangle = triangle + ntr;
    ptriangle->A = A; ptriangle->B = B; ptriangle->C = C;
    ptriangle->a = a; ptriangle->b = b; ptriangle->c = c;
    ptriangle->h = h;

    XA = xA / zA; YA = yA / zA;
    XB = xB / zB; YB = yB / zB;
    XC = xC / zC; YC = yC / zC;

    DA = XB * YC - XC * YB;
    DB = XC * YA - XA * YC;
    DC = XA * YB - XB * YA;
    D = DA + DB + DC;
    DAB = DC - M * (XA - XB); DAC = DB - M * (XC - XA); DBC = DA - M * (XB - XC);

    topcode[0] = (D * DAB > 0); topcode[1] = (D * DAC  > 0); topcode[2] = (D * DBC > 0);
    Xleft[0] = XA; Yleft[0] = YA; Xright[0] = XB; Yright[0] = YB;
    Xleft[1] = XA; Yleft[1] = YA; Xright[1] = XC; Yright[1] = YC;
    Xleft[2] = XB; Yleft[2] = YB; Xright[2] = XC; Yright[2] = YC;

    int i;
    for (i = 0; i < 3; ++i)
        if (Xleft[i] > Xright[i] || (Xleft[i] == Xright[i] && Yleft[i] > Yright[i])) {
            aux = Xleft[i]; Xleft[i] = Xright[i]; Xright[i] = aux;
            aux = Yleft[i]; Yleft[i] = Yright[i]; Yright[i] = aux;
        }
    ipixmin = xwhole(min3(XA, XB, XC));
    ipixmax = xwhole(max3(XA, XB, XC));

    for (ipix = ipixmin; ipix <= ipixmax; ++ipix) {
        lower[ipix] = up[ipix] = 10000;
        upper[ipix] = low[ipix] = -10000;
    }

    for (i = 0; i < 3; ++i) {
        ipixleft = xwhole(Xleft[i]); ipixright = xwhole(Xright[i]);
	denom = Xright[i] - Xleft[i];
        if (ipixleft != ipixright)
            slope = (Yright[i] - Yleft[i]) / denom;

        j_old = ywhole(Yleft[i]);
        for (ipix = ipixleft; ipix <= ipixright; ++ipix) {
            if (ipix == ipixright)
                jl = ywhole(Yright[i]);
            else
                jl = ywhole(Yleft[i] + (xreal(ipix + 1) - Xleft[i]) * slope);
            if (topcode[i]) {
                upper[ipix] = max3(j_old, jl, upper[ipix]);
                up[ipix] = min3(j_old, jl, up[ipix]);
            } else {
                lower[ipix] = max3(j_old, jl, lower[ipix]);
                low[ipix] = min3(j_old, jl, low[ipix]);
            }
            j_old = jl;
        }
    }

    for (ipix = ipixmin; ipix <= ipixmax; ++ipix)
        for (jpix = lower[ipix]; jpix <= upper[ipix]; ++jpix) {
            pointer = &(ScreeN[ipix][jpix]);
            if (jpix > low[ipix] && jpix < up[ipix]) {
                xR = Xmin + (ipix + 0.5) * deltaX;
                yR = Ymin + (jpix + 0.5) * deltaY;
                denom = a * xR + b * yR + c * d;
                dist = fabs(denom) > eps ? h * sqrt(xR*xR + yR*yR + 1.) / denom : BIG;

                if (dist < pointer->tr_dist) {
                    pointer->tr_cov = ntr;
                    pointer->tr_dist = dist;
                }
            } else {
                pnode = (struct node *)malloc(sizeof(struct node));
                if (pnode == NULL)
                    error("Allocation error \n");
                pnode->jtr = ntr; pnode->next = pointer->start;
                pointer->start = pnode;
            }
        }
    ntr++;
}

