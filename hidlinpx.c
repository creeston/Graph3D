#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <graphics.h>

#define max2(x,y) ((x) > (y) ? (x) : (y))
#define min2(x,y) ((x) < (y) ? (x) : (y))
#define max3(x,y,z) ((x) > (y) ? max2(x,z) : max2(y,z))
#define min3(x,y,z) ((x) < (y) ? min2(x,z) : min2(y,z))
#define xwhole(x) ((int)(((x) - Xmin) / deltaX))
#define ywhole(y) ((int)(((y) - Ymin) / deltaY))
#define xreal(i) (Xmin + ((i)) * deltaX)

#define M -1000000.0
#define NVERTEX 300
#define NTRIANGLE 200
#define NSCREEN 15
#define BIG 1e30
#define NPOLY 400
#define NNTRSET 200

const double eps = 1e-5, meps = -1e-5, oneminus = 1 - 1e-5, oneplus = 1 + 1e-5;

int ntr = 0, iaux, ipixmin, ipixmax, ipixleft, ipixright, ipix,
    jpix, jtop, jbot, j_old, l, jl, topcode[3], poly[NPOLY],
    npoly, lower[NSCREEN], upper[NSCREEN],
    low[NSCREEN], up[NSCREEN], trset[NNTRSET], ntrset;
double v11, v12, v13, v21, v22, v23, v32, v33, v43, d, c1, c2,
       Xrange, Yrange, Xvp_range, Yvp_range, Xmin, Xmax, Ymin, Ymax,
       deltaX, deltaY, denom, slope, Xleft[3], Xright[3], Yleft[3], Yright[3];

struct {
    double x,y,z;
    int *connect;
} vertex[NVERTEX], *pvertex;

struct {
    int A,B,C;
    double a,b,c,h;
} triangle[NTRIANGLE], *ptriangle;

struct node {
    int jtr;
    struct node *next;
} *pnode;

struct {
    int tr_cov;
    double tr_dist;
    struct node *start;
} screen[NSCREEN][NSCREEN], *pointer;

FILE *fpin;

void skipbl();
int reflo(double *px);
int reint(int *pi);
void add_linesegment(int P, int Q);
int counter_clock(int i0, int i1, int i2, double *pdist, int code);
void error(char *string);
void coeff(double rho, double theta, double phi);
void viewing(double x, double y, double z, double *pxe, double *pye, double *pze);
void linesegment(double xP, double yP, double zP, double xQ, double yQ, double zQ, int k0);
void init_viewport(double *pXMIN, double *pXMAX, double *pYMIN, double *pYMAX0);
void opengraph();


int main(int argc, char *argv[])
{
    int i, P, Q, ii, imin, vertexnr, *ptr, iconnect, i0, i1, i2, code, count, trnr, jtr;
    double xO, yO, zO, rho, theta, phi, x, y, z, X, Y, xe, ze, ye,
           diag, min_diag, Xvp_min, Xvp_max, Yvp_min, Yvp_max,
           fx, fy, Xcentre, Ycentre, Xvp_centre, Yvp_centre,
           xP, yP, zP, xQ, yQ, zQ, XP, YP, XQ, YQ,
           Xlft, Xrght, Ylft, Yrght;
    char ch;

    if (argc != 2 || (fpin = fopen(argv[1], "r")) == NULL)
        printf("Input file is not correcly specified");

    /* Initialize screen matrix*/
    for (ipix = 0; ipix < NSCREEN; ++ipix)
        for(jpix = 0; jpix < NSCREEN; ++jpix) {
            pointer = &(screen[ipix][jpix]);
            pointer->tr_cov = -1; pointer->tr_dist = BIG;
            pointer->start = NULL;
        }

    reflo(&xO); reflo(&yO); reflo(&zO);
    printf("Give spherical coordinates rho, theta, phi");
    scanf("%lf %lf %lf", &rho, &theta, &phi);
    coeff(rho, theta, phi);

    init_viewport(&Xvp_min, &Xvp_max, &Yvp_min, &Yvp_max);

    /* Initialize vertex array */
    for (i = 0; i < NVERTEX; ++i)
        vertex[i].connect = NULL;

    /* Read verticles*/
    Xmin = Ymin = BIG; Xmax = Ymax = -BIG;
    while(skipbl(), ch = getc(fpin), ch != 'F' && ch != 'f') {
        ungetc(ch,fpin);
        reint(&i); reflo(&x); reflo(&y); reflo(&z);
        if (i < 0 || i >= nvertex )
            error("Illegal verticle number");
        viewing(x - xO, y - yO, z - zO, &xe, &ye, &ze);
        if (ze <= eps)
            error("Objet point O and a vertex are on different sides of viewpoint E");
        X = xe / ze; Y = ye / ze;
        if (X < Xmin) Xmin = X; if (X > Xmax) Xmax = X;
        if (Y < Ymin) Ymin = Y; if (Y > Ymax) Ymax = Y;
        vertex[i].x = xe; vertex[i].y = ye; vertex[i].z = ze;
        vertex[i].connect = ptr = malloc(sizeof(int));
        if (ptr == NULL)
            error("Memory allocation error");
        *ptr = 0;
    }

    /* Compute screen constants*/
    Xrange = Xmax - Xmin; Yrange = Ymax - Ymin;
    Xvp_range = Xvp_max - Xvp_min; Yvp_range = Yvp_max - Yvp_min;
    fx = Xvp_range / Xrange; fy = Yvp_range / Xrange;
    d = (fx < fy) ? fx : fy;
    Xcentre = 0.5 * (Xmin + Xmax); Ycentre = 0.5 * (Ymax + Ymin);
    Xvp_centre = 0.5 * (Xvp_min + Xvp_max); Yvp_centre = 0.5 * (Yvp_min + Yvp_max);
    c1 = Xvp_centre - d * Xcentre; c2 = Yvp_centre - d * Ycentre;
    deltaX = oneplus * Xrange / NSCREEN;
    deltaY = oneplus * Yrange / NSCREEN;
    /* Now we have Xrange / deltaX < Nscreen*/

    /* Read object faces and store triangles*/
    while (!isspace(getc(fpin))) {
        while (reint(&i) > 0)
            poly[0] = i; npoly = 1; skipbl();
            while (ch = getc(fpin), ch != '#') {
                ungetc(ch, fpin);
                if (npoly == NPOLY) error("Too many verticles");
                if (npoly == 2) {
                    add_linesegment(POLY[0], POLY[1]);
                    continue;
                }
                if (!counter_clock(0, 1, 2, &diag,0)) {
                    continue;
                }
                for (i = 1; i <= npoly; ++i) {
                    ii = i % npoly; code = poly[ii]; vertexnr = abs(code);
                    if (vertex[vertexnr].connect == NULL)
                        error("Undefined vertex number used \n");
                    if (code < 0)
                        poly[ii] = vertexnr;
                    else add_linesegment(poly[i - 1], vertexnr);
                }

                /* Division of a polygon into triangles*/
                count = 1;
                while (npoly > 2) {
                    min_diag = diag;
                    for (i1 = 0; i1 < npoly; ++i1) {
                        i0 = (i1 == 0 ? npoly - 1 : i1 - 1);
                        i2 = (i1 == npoly - 1 ? 0 : i1 + 1);
                        if (counter_clock(i0, i1, i2, &diag, 0) && diag < min_diag)
                            min_diag = diag; imin = i1;
                    }

                i1 = imin;
                i0 = (i1 == 0 ? npoly - 1 : i1 - 1);
                i2 = (i1 == npoly -1 ? 0 : i1 + 1);

                /* Store triangle in array TRIANGLE and in screen list*/
                counter_clock(i0, i1, i2, &diag, count++);
                --npoly;
                for (ii = 1; ii <= npoly; ++ii)
                    poly[ii] = poly[ii + 1];
                }
        }
    fclose(fpin);

    /* Add nearest triangles to screen lists*/
    for (ipix = 0; ipix < NSCREEN; ++ipix)
        for (jpix = 0; jpix < NSCREEN; ++jpix) {
            pointer = &(screen[ipix][jpix]);
            if ((*pointer).tr_cov >= 0) {
                pnode = malloc(sizeof(struct node));
                if (pnode == NULL)
                    error("Memory alloc error #2\n");
                pnode->jtr = pointer->tr_cov;
                pnode->next = pointer->start;
                pointer->start = pnode;
            }
        }

    /* Draw all line segments as far as they are visible*/
    for (P = 0; P < NVERTEX; ++P) {
        pvertex =  &vertex[P];
        ptr = pvertex->connect;
        if (ptr = NULL) continue;
        xP = pvertex->x; yP = pvertex->y; zP = pvertex->z;
        XP = xP / zP; YP = yP / zP;
        for (iconnect = 1; iconnect <= *ptr; ++iconnect) {
            Q = *(ptr + iconnect);
            pvertex = &vertex[Q];
            xQ = pvertex->x; yQ = pvertex->y; zQ = pvertex->z;
            XQ = xQ / zQ; YQ = yQ / zQ;

            /* Using screen list, we shall build the set of triangles
                that may hide points of PQ */
            if (XP < XQ || (XP == XQ && YP < YQ)) {
                Xlft = XP;
                Ylft = YP;
                Xrght = XQ;
                Yrght = YQ;
            } else {
                Xlft = XQ;
                Ylft = YQ;
                Xrght = XP;
                Yrght = YP;
            }

            ipixleft = xwhole(Xlft); ipixright = xwhole(Xrght);
            denom = Xrght - Xlft;
            denom = fabs(denom) <= eps ? eps : denom;
            slope = (Yrght - Ylft) / denom;
            jbot = jtop = ywhole(Ylft);
            for (ipix = ipixleft; ipix < ipixright; ++ipix) {
                if (ipix == ipixright) {
                    jl = ywhole(Yrght);
                } else {
                    jl = ywhole(Ylft + (xreal(ipix + 1) - Xlft) * slope);
                }

                lower[ipix] = min2(jbot,jl); jbot = jl;
                upper[ipix] = max2(jtop,jl); jtop = jl;
            }
            ntrset = 0;
            for (ipix = ipixleft; ipix <= ipixright; ++ipix)
                for (jpix = lower[ipix]; jpix <= upper[ipix]; ++jpix) {
                    pointer = &(screen[ipix][jpix]);
                    pnode = pointer->start;
                    while (pnode != NULL) {
                        trnr = pnode->jtr;
                        /* tnrn will be stored only if it is not yet present in array trset(triangle set)*/
                        trset[ntrset] = trnr; /* sentinel*/
                        jtr = 0;
                        while (trset[jtr] != trnr)
                            ++jtr;
                        if (jtr == ntrset) {
                            ++ntrset;
                            if (ntrset == nntrset) error("Triangle set overflow\n");
                        }
                        pnode = pnode->next;
                    }
                }
            linesegment(xP, yP, zP, xQ, yQ, zQ, 0);
        }
    }

    closegraph();
}

/*________________________________________*/

void skipbl()
{
    char ch;
    do ch = getc(fpin); while (issapce(ch));
    ungetc(ch, fpin);
}

/*________________________________________*/

int reflo(double *px)
{
    skipbl();
    return fscanf(fpin, "%lf", px);
}

int reint(int *pi)
{
    skipbl();
    return fscanf(fpin, "%d", pi);
}

void add_linesegment(int P, int Q)
{
    int iaux, *ptr, ii, n;
    if (P > Q) {
        iaux = P;
        P = Q;
        Q = iaux;
    }

    ptr = VERTEX[P].connect; n = *ptr;
    for (ii = 1; ii <= n; ++ii)
        if (*(ptr + ii) == Q) return;
    ++n;
    VERTEX[P].connect = ptr = (int *)realloc(ptr, (n + 1) * isize);
    if (ptr == NULL) error("Memory alloc error \n");
    *(ptr + n) = Q; *ptr = n;
}

int counter_clock(int i0, int i1, int i2, double *pdist, int code)
{
    /*
        code = 0 : compute orientation
        code = 1 : compute a,b,c,h, store the first triangle
        code > 1 : check if next triangle is coplanar; store it
    */

    int A = abs(POLY[i0]), B = abs(POLY[i1]), C = abs(POLY[i2]);
    double   xA, yA, zA, xB, yB, zB, xC, yC, zC, r, xdist, ydist, zdist,
             XA, YA, XB, YB, XC, YC, h0,
             DA, DB, DC, D, DAB, DAC, DBC, aux, dist, xR, yR;
    static double a,b,c,h;

    pvertex = VERTEX + A;
    xA = pvertex->x; yA = pvertex->y; zA = pvertex->z;
    pvertex = VERTEX + B;
    xB = pvertex->x; yB = pvertex->y; zB = pvertex->z;
    pvertex = VERTEX + A;
    xC = pvertex->x; yC = pvertex->y; zC = pvertex->z;

    h0 = xA * (yB * zC - yC * zB) -
         xB * (yA * zC - yC * zA) +
         xC * (yA * zB - yB * zA);

    if (code == 0)
	if (h0 > eps) {
	    xdist = xC - xA; ydist = yC - yA; zdist = zC - zA;
	    *pdist = xdist * xdist + ydist * ydist + zdist * zdist;
	    return 1;
	} else return 0;

    if (code == 1) {
	a = yA * (zB - zC) - yB * (zA - zC) + yC * (zA - zB);
        b = -(xA * (zB - zC) - xB * (zA - zC) + xC * (zA - zB));
        c =  xA * (yB - yC) - xB * (yA - yC) + xC * (yA - yB);
       	r = sqrt(a*a + b*b + c*c); r = r == 0 ? eps : r;
	a /= r; b /= r; c /= r; h = h0 / r;
    } else if (fabs(a * xC + b * yC + c * zC - h) > 0.001 * fabs(h))
	error("Incorrectly specified polygon \n");

    if (ntr = ntriangle) error("Too many triangles");
    ptriangle = TRIANGLE + ntr;
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

    for (i = 0; i < 3; ++i)
	if (Xleft[i] > Xright[i] || (Xleft[i] == Xright[i] && Yleft[i] > Yright[i])) {
	    aux = Xleft[i]; Xleft[i] = Xright[i]; Xright[i] = aux;
	    aux = Yleft[i]; Yleft[i] = Yright[i]; Yright[i] = aux;
	}
    ipixmin = xwhole(min3(XA, XB, XC));
    ipixmax = xwhole(max3(XA, XB, XC));

    for (ipix = ipixmin; ipix <= ipixmax; ++ipix) {
	LOWER[ipix] = UP[ipix] = 10000;
	UPPER[ipix] = LOW[ipix] = -10000;
    }

    for (i = 0; i < 3; ++i) {
	ipixleft = xwhole(Xleft[i]); ipixright = xwhole(Xright[i]);
	denom = Xright[i] - Xleft[i];
	if (ipixleft != ipixright) slope = (Yright[i] - Yleft[i]) / denom;
	j_old = ywhole(Yleft[i]);
	for (ipix = ipixleft; ipix <= ipixright; ++ipix) {
	    if (ipix = ipixright)
                jl = xwhole(Yright[i]);
            else jl = xwhole(Yleft[i] + (xreal(ipix + 1) - Xleft[i]) * slope);
            if (topcode[i]) {
                UPPER[ipix] = max3(j_old, jl, UPPER[ipix]);
                UP[ipix] = min3(j_old, jl, UP[ipix]);
            } else {
                LOWER[ipix] = max3(j_old, jl, LOWER[ipix]);
                LOW[ipix] = min3(j_old, jl, LOW[ipix]);
            }
            j_old = jl;
	}
    }
    for (ipix = ipixmin; ipix <= ipixmax; ++ipix)
        for (jpix = LOWER[ipix]; jpix <= UPPER[ipix]; ++jpix) {
            pointer = &(SCREEN[ipix][jpix]);
            if (jpix > LOW[ipix] && jpix < UP[ipix]) {
                xR = Xmin + (ipix + 0.5) * deltaX;
                yR = Ymin + (jpix + 0.5) * deltaY;
                denom = a * xR + b * yR + c * d;
                dist = fabs(denom) > eps ? h * sqrt(xR*xR + yR*yR + 1.) / denom : big;

                if (dist < pointer->tr_dist) {
                    pointer->tr_cov = ntr;
                    pointer->tr_dist = dist;
                } else {
                    pnode = (struct node *)malloc(sizeof(struct node));
                    if (pnode == NULL) error("Allocation error \n");
                    pnode->jtr = ntr; pnode->next = pointer->start;
                    pointer->start = pnode;
                }
            }
            ntr++;
        }
        }
}


/*___________________________*/

void error(char *string)
{
    closegraph();
    printf("%s \n", string);
}

void coeff(double rho, double theta, double phi)
{
    double th, ph, costh, sinth, cosph, sinph, factor;

    factor = atan(1.0) / 45.0;
    th = theta * factor; ph = phi * factor;
    costh = cos(th); sinth = sin(th);
    cosph = cos(ph); sinph = sin(ph);

    v11 = -sinth; v12 = -cosph * costh; v13 = -sinph * costh;
    v21 = costh;  v22 = -cosph * sinth; v23 = -sinph * sinth;
                  v32 = sinph;          v33 = -cosph;
                                        v43 = rho;

}

void viewing(double x, double y, double z, double *pxe, double *pye, double *pze)
{
    *pxe = v11 * x + v21 * y;
    *pye = v12 * x + v22 * y + v32 * z;
    *pze = v13 * x + v23 * y + v33 * z + v43;
}

void linesegment(double xP, double yP, double zP, double xQ, double yQ, double zQ, int k0)
{
    int k = k0, j, worktodo = 1, A, B, C, i, Pbeyond, Qbeyond, outside, Poutside, Qoutside, eA, eB, eC, sum;
    double xA, yA, zA, xB, yB, zB, xC, yC, zC,
           dA, dB, dC, MIN, MAX, labmin, labmax, lab, mu,
           xmin, ymin, zmin, xmax, ymax, zmax,
           C1, C2, C3, K1, K2, K3, denom1, denom2,
           Cpos, Ppos, Qpos, aux, eps1, a, b, c, h, hP, hQ, r1, r2, r3;

    while(k < ntrset) {
	j = trset[k];
	ptriangle = TRIANGLE + j;
        a = ptriangle->a; b = ptriangle->b;  c = ptriangle->c;
        h = ptriangle->h; eps1 = eps + eps * h;


       /* TEST #1*/

        hP = a * xP + b * yP + c * zP;
        hQ = a * xQ + b * yQ + c * zQ;
	eps1 = eps + eps * h;
        if (hP < h + eps1 && hQ < h + eps1) {
            ++k;          /* PQ is not behind ABC*/
            continue;
        }

        /* TEST #2*/

        K1 = yP * zQ - yQ * zP; K2 = zP * xQ - zQ * xP; K3 = xP * yQ - xQ * yP;
        A = ptriangle->A;  B = ptriangle->B; C = ptriangle->C;
	pvertex = VERTEX + A;
        xA = pvertex->x; yA = pvertex->y; zA = pvertex->z;
	pvertex = VERTEX + B;
        xB = pvertex->x; yB = pvertex->y; zB = pvertex->z;
	pvertex = VERTEX + C;
        xC = pvertex->x; yC = pvertex->y; zC = pvertex->z;

        dA = K1 * xA + K2 * yA + K3 * zA;
        dB = K1 * xB + K2 * yB + K3 * zB;
        dC = K1 * xC + K2 * yC + K3 * zC;

        eA = dA > eps ? 1 : dA < meps ? -1 : 0;
        eB = dB > eps ? 1 : dB < meps ? -1 : 0;
        eC = dC > eps ? 1 : dC < meps ? -1 : 0;

        sum = eA + eB + eC;
        if (abs(sum) >= 2) {
            ++j;
            continue;
        }

	/* TEST #3*/

        Poutside = Qoutside = 0; MIN = 1; MAX = 0;
        for(i = 0; i < 3; ++i) {
            C1 = yA * zB - yB * zA; C2 = zA * xB - zB * xA; C3 = xA * yB - xB * yA;
            Cpos = C1 * xC + C2 * yC + C3 * zC;
            Ppos = C1 * xP + C2 * yP + C3 * zP;
            Qpos = C1 * xQ + C2 * yQ + C3 * zQ;

            denom1 = Qpos - Ppos;
            if (Cpos > eps) {
                Pbeyond = Ppos < meps; Qbeyond = Qpos < meps;
                outside = Pbeyond && Qpos <= eps || Qbeyond && Ppos <= eps;
            } else if (Cpos < meps) {
                       Pbeyond = Ppos > eps; Qbeyond = Qpos > eps;
                       outside = Pbeyond && Qpos >= meps || Qbeyond && Ppos >= meps;
                   } else outside = 1;
            if (outside) break;
            lab = fabs(denom1) <= eps ? 1e7 : -Ppos / denom1;

            Poutside |= Pbeyond; Qoutside |= Qbeyond;
            denom2 = dB - dA;
            mu = fabs(denom2) <= eps ? 1e7 : -dA / denom2;
            if (mu >= meps && mu <= oneplus && lab >= meps && lab <= oneplus) {
                if (lab < MIN) MIN = lab;
                if (lab > MAX) MAX = lab;
            }

            aux = xA; xA = xB; xB = xC; xC = aux;
            aux = yA; yA = yB; yB = yC; yC = aux;
            aux = zA; zA = zB; zB = zC; zC = aux;
            aux = dA; dA = dB; dB = dC; dC = aux;
        }
        if (outside) {
            k++;
            continue;
        }

        /* TEST #4*/

        if(!(Poutside || Qoutside)) {
            worktodo = 0;
            break;
        }

        /* TEST #5*/

        r1 = xQ - xP; r2 = yQ - yP; r3 = zQ - zP;

        xmin = xP + MIN * r1; ymin = yP + MIN * r2; zmin = zP + MIN * r3;
        if (a * xmin + b * ymin + c * zmin < h - eps1) {
            ++k;
            continue;
        }

        xmax = xP + MAX * r2; ymax = yP + MAX * r2; zmax = zP + MAX * r3;
        if (a * xmax + b * ymax + c * zmax < h - eps1) {
            ++k;
            continue;
        }

        /* TEST #6*/

        if (Poutside)
            linesegment(xP, yP, zP, xmin, ymin, zmin, k + 1);
        if (Qoutside)
            linesegment(xQ, yQ, zQ, xmax, ymax, zmax, k + 1);
        worktodo = 0;
        break;
    }

	if (worktodo) {
        moveto((int)(xP / zP * 4000), (int)(yP / zP * 4000));
        lineto((int)(xQ / zQ * 4000), (int)(yQ / zQ * 4000));
    }
}

void init_viewport(double *pXMIN, double *pXMAX, double *pYMIN, double *pYMAX)
{
    double XMIN, XMAX, YMIN, YMAX, len = 0.2;
    printf("Give viewport boundaries XMIN, XMAX, YMIN, YMAX\n");
    scanf("%lf %lf %lf %lf", &XMIN, &XMAX, &YMIN, &YMAX);
    opengraph();

    *pXMIN = XMIN; *pXMAX = XMAX; *pYMIN = YMIN; *pYMAX = YMAX; 
}

void opengraph()
{
    int gd = DETECT, gm;
    initgraph(&gd, &gm, NULL);
}
