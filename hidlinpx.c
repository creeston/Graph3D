#include<stdio.h>
#include<math.h>
#include<ctype.h>
#include<graphics.h>
#include "test_header.h"
#include<X11/Xlib.h>
#include<time.h>
#include<SDL/SDL.h>


int ntr, iaux, ipixmin, ipixmax, ipixleft, ipixright, ipix,
    jpix, jtop, jbot, j_old, l, jl, topcode[3], poly[NPOLY],
    npoly, lower[NSCREEN], upper[NSCREEN],
    low[NSCREEN], up[NSCREEN], trset[NNTRSET], ntrset, isize = sizeof(int);
double v11, v12, v13, v21, v22, v23, v32, v33, v43, d, c1, c2,
       Xrange, Yrange, Xvp_range, Yvp_range, Xmin, Xmax, Ymin, Ymax,
       deltaX, deltaY, denom, slope, Xleft[3], Xright[3], Yleft[3], Yright[3];

/*struct ver {
    double x,y,z;
    int *connect;
};

struct ver vertex[NVERTEX], *pvertex;

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
} Screen[NSCREEN][NSCREEN], *pointer;
*/

struct ver vertex[NVERTEX], *pvertex;
struct tr triangle[NTRIANGLE], *ptriangle;
struct node *pnode;
struct scr ScreeN[NSCREEN][NSCREEN], *pointer;

FILE *fpin;

void drawsegment(double xP, double yP, double zP, double xQ, double yQ, double zQ);
void polydiv();
void update();
void check();
void update();
void border_draw(double XMIN, double XMAX, double YMIN, double YMAX);
void drawtr();
void extra_init(double rho, double theta, double phi, double xMIN, double xMAX, double yMIN, double yMAX);
void read_vertex();
void scale_calc(double Xvp_min, double Xvp_max, double Yvp_min, double Yvp_max);
void read_faces();
void tr_add();
void line_and_pix();
void read_and_draw();


int main(int argc, char *argv[])
{
    int i, running = 1, rotate = 0, xcountpos = 0, xcountneg = 0, ycountpos = 0, ycountneg = 0;
    double Xvp_min, Xvp_max, Yvp_min, Yvp_max, rho = 3000, theta = 45, phi = 0;
    char ch;
    int dx, dy, x, y;


    init_viewport(&Xvp_min, &Xvp_max, &Yvp_min, &Yvp_max);
    opengraph();

    /* Открытие файла*/
    if (argc != 2 || (fpin = fopen(argv[1], "r")) == NULL)
        printf("Input file is not correcly specified");
    ntr = 0;

    read_and_draw(rho, theta, phi, Xvp_min, Xvp_max, Yvp_min, Yvp_max);

    while (running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
	    switch (event.type) {
		case SDL_QUIT:
                    running = 0;
		    closegraph();
		    return 0;
                    break;
                case SDL_MOUSEBUTTONDOWN:
                    x = event.button.x;
                    y = event.button.y;
                    if (x > Xvp_min && x < Xvp_max && y > Yvp_min && y < Yvp_max / 2)
                        rotate = 1;
                    if (x > Xvp_min && x < Xvp_max / 2 && y > Yvp_min && y < Yvp_max)
                        rotate = 1;
                    break;
                case SDL_MOUSEBUTTONUP:
                    rotate = 0;
                    break;
                case SDL_MOUSEMOTION:
                    if (rotate && event.motion.x < Xvp_max && event.motion.y < Yvp_max) {
                        dx = event.motion.xrel;
                        dy = event.motion.yrel;
                        if (dy > 0) ycountpos++;
                        if (dy < 0) ycountneg++;
                        if (dx > 0) xcountpos++;
                        if (dx < 0) xcountneg++;

                        if (ycountpos == 4) {
                            phi += 2;
                            update();
                            read_and_draw(rho, theta, phi, Xvp_min, Xvp_max, Yvp_min, Yvp_max);
                            ycountpos = 0;
                        }
                        if (ycountneg == 4) {
                            phi -= 2;
                            update();
                            read_and_draw(rho, theta, phi, Xvp_min, Xvp_max, Yvp_min, Yvp_max);
                            ycountneg = 0;
                        }
                        if (xcountpos == 4) {
                            theta +=2;
                            update();
                            read_and_draw(rho, theta, phi, Xvp_min, Xvp_max, Yvp_min, Yvp_max);
                            xcountpos = 0;
                        }
                        if (xcountneg == 4) {
                            theta -= 2;
                            update();
                            read_and_draw(rho, theta, phi, Xvp_min, Xvp_max, Yvp_min, Yvp_max);
                            xcountneg = 0;
                        }

                        //update();
                        //read_and_draw(rho, theta, phi, Xvp_min, Xvp_max, Yvp_min, Yvp_max);

                        break;
                    }
             }

        }
    }
}

void check()
{
   printf("checkpoint\n");
   char c;
}

void update()
{
    int i;

    delay(20);
    cleardevice();

    for (i = 0; vertex[i].connect != NULL; ++i)
        free(vertex[i].connect);

    for (ipix = 0; ipix < NSCREEN; ++ipix)
        for (jpix = 0; jpix < NSCREEN; ++jpix)
            free(ScreeN[ipix][jpix].start);
    ntr = 0;
    fseek(fpin, 0L, SEEK_SET);

}

void drawsegment(double xP, double yP, double zP, double xQ, double yQ, double zQ)
{
    moveto((int)(d * xP / zP + c1), (int)(yP * d / zP + c2));
    lineto((int)(xQ * d / zQ + c1), (int)(yQ * d / zQ + c2));
}

void drawtr()
{
    int i, XA, YA, XB, YB, XC, YC, col;
    for (i = 0; i < ntr; ++i) {
        col = (int)rand() % 10 + 1;
        setcolor(RED);

        XA = (int)(vertex[triangle[i].A].x * d / vertex[triangle[i].A].z + c1);
	YA = (int)(vertex[triangle[i].A].y * d / vertex[triangle[i].A].z + c2);
	XB = (int)(vertex[triangle[i].B].x * d / vertex[triangle[i].B].z + c1);
	YB = (int)(vertex[triangle[i].B].y * d / vertex[triangle[i].B].z + c2);
	XC = (int)(vertex[triangle[i].C].x * d / vertex[triangle[i].C].z + c1);
	YC = (int)(vertex[triangle[i].C].y * d / vertex[triangle[i].C].z + c2);
        moveto(XA, YA);
	lineto(XB, YB);
        lineto(XC, YC);
	lineto(XA, YA);
	floodfill((XA + XB + 2 * XC) / 4, (YA + YB + 2 * YC) / 4, RED);
    }
    update();

}

void border_draw(double XMIN, double XMAX, double YMIN, double YMAX)
{
    int len = 5;

    moveto(XMIN, YMIN + len); lineto(XMIN, YMIN); lineto(XMIN + len, YMIN);
    moveto(XMAX - len, YMIN); lineto(XMAX, YMIN); lineto(XMAX, YMIN + len);
    moveto(XMAX, YMAX - len); lineto(XMAX, YMAX); lineto(XMAX - len, YMAX);
    moveto(XMIN + len, YMAX); lineto(XMIN, YMAX); lineto(XMIN, YMAX - len);
    moveto((XMIN + XMAX) / 2, YMIN); lineto((XMIN + XMAX) / 2, YMIN);
}


void polydiv()
{
    int count, i1, i2, i0, ii, imin;
    double diag, min_diag;
    count = 1;
    while (npoly > 2) {
        min_diag = BIG;
        for (i1 = 0; i1 < npoly; ++i1) {
            i0 = (i1 == 0 ? npoly - 1 : i1 - 1);
            i2 = (i1 == npoly - 1 ? 0 : i1 + 1);
            if (counter_clock(i0, i1, i2, &diag, 0) && diag < min_diag) {
                min_diag = diag;
                imin = i1;
            }
        }

        i1 = imin;
        i0 = (i1 == 0 ? npoly - 1 : i1 - 1);
        i2 = (i1 == npoly - 1 ? 0 : i1 + 1);

        /* Запись треугольника в массив triangle и ассоциация его с пикслями */
        counter_clock(i0, i1, i2, &diag, count++);
        --npoly;
        for (ii = i1; ii <= npoly; ++ii)
            poly[ii] = poly[ii + 1];
     }

}

void extra_init(double rho, double theta, double phi, double xMIN, double xMAX, double yMIN, double yMAX)
{
    int i;

    /* Инициализируем экранную матрицу*/
    for (ipix = 0; ipix < NSCREEN; ++ipix)
        for(jpix = 0; jpix < NSCREEN; ++jpix) {
            pointer = &(ScreeN[ipix][jpix]);
            pointer->tr_cov = -1; pointer->tr_dist = BIG;
            pointer->start = NULL;
        }

    /* Считаем матрицу поворота*/
    coeff(rho, theta, phi);

    /* задаем область вывода */
    border_draw(xMIN, xMAX, yMIN, yMAX);

    /* Инициализируем массив вершин */
    for (i = 0; i < NVERTEX; ++i)
        vertex[i].connect = NULL;


}

void read_vertex()
{
    int i, *ptr;
    char ch;
    double x, y, z, xO, yO, zO, xe, ye, ze, X, Y;

    fscanf(fpin, "%lf %lf %lf", &xO, &yO, &zO);
    Xmin = Ymin = BIG; Xmax = Ymax = -BIG;
    while(skipbl(), ch = getc(fpin), ch != 'F' && ch != 'f') {
        ungetc(ch,fpin);
        //reint(&i); reflo(&x); reflo(&y); reflo(&z);
        fscanf(fpin, "%d %lf %lf %lf", &i, &x, &y, &z);
        if (i < 0 || i >= NVERTEX )
            error("Неверный номер вершины \n");
        //преобразуем координаты (через матрицу поворота)
        viewing(x - xO, y - yO, z - zO, &xe, &ye, &ze);
        if (ze <= eps)
            error("");

        //находим max и min для масштабирования в наше окно вывода
        X = xe / ze; Y = ye / ze;
        if (X < Xmin) Xmin = X; if (X > Xmax) Xmax = X;
        if (Y < Ymin) Ymin = Y; if (Y > Ymax) Ymax = Y;

        vertex[i].x = xe; vertex[i].y = ye; vertex[i].z = ze;
        vertex[i].connect = ptr = (int *)malloc(sizeof(int));
        if (ptr == NULL)
            error("Memory allocation error");
        *ptr = 0;
    }

}

void scale_calc(double Xvp_min, double Xvp_max, double Yvp_min, double Yvp_max)
{
    double fx, fy, Xcentre, Ycentre, Xvp_centre, Yvp_centre;

    Xrange = Xmax - Xmin; Yrange = Ymax - Ymin;
    Xvp_range = Xvp_max - Xvp_min; Yvp_range = Yvp_max - Yvp_min;
    fx = Xvp_range / Xrange; fy = Yvp_range / Yrange;
    d = (fx < fy) ? fx : fy;
    Xcentre = (Xmin + Xmax) / 2; Ycentre = (Ymax + Ymin) / 2;
    Xvp_centre = (Xvp_min + Xvp_max) / 2; Yvp_centre = (Yvp_min + Yvp_max) / 2;
    c1 = Xvp_centre - d * Xcentre; c2 = Yvp_centre - d * Ycentre;
    deltaX = oneplus * Xrange / NSCREEN; // = Xrange / Nscreen
    deltaY = oneplus * Yrange / NSCREEN;

}

void read_faces()
{
    int i, ii;
    char ch;
    double diag;
    while (!isspace(getc(fpin)))
        ;

    while (reint(&i) > 0) {
        poly[0] = i; npoly = 1; skipbl();
        while (ch = getc(fpin), ch != '#') {
            ungetc(ch, fpin);
            reint(&poly[npoly++]);
            if (npoly == NPOLY)
                error("Слишком много вершин\n");
        }

        if (npoly == 1)
            error("Только одна вершина\n");

        if (npoly == 2) {
           add_linesegment(poly[0], poly[1]);
           continue;
        }

        //проверяем, является ли треугольник задним
        if (!(counter_clock(0, 1, 2, &diag, 0)))
            continue;

        //записываем в поле connect информацию об соединениях вершин -> формируем отрезки PQ
        for (i = 1; i <= npoly; ++i) {
            ii = i % npoly;

            if (vertex[abs(poly[ii])].connect == NULL)
               error("Неопределенная вершина \n");

            if (poly[ii] < 0)
                poly[ii] = abs(poly[ii]);
            else add_linesegment(poly[i - 1], poly[ii]);
        }

        /* Разбиение треугольника на полигоны */
        polydiv();
    }
//!//    fclose(fpin);
}

void tr_add()
{
    for (ipix = 0; ipix < NSCREEN; ++ipix)
        for (jpix = 0; jpix < NSCREEN; ++jpix) {
            pointer = &(ScreeN[ipix][jpix]);

            if ((*pointer).tr_cov >= 0) {
                pnode = (struct node *)malloc(sizeof(struct node));
                if (pnode == NULL)
                    error("Memory alloc error #2\n");

                pnode->jtr = pointer->tr_cov;
                pnode->next = pointer->start;
                pointer->start = pnode;
            }
        }

}

void line_and_pix()
{
    int P,Q, iconnect, *ptr, i, jtr, trnr;
    double xQ, yQ, zQ, xP, yP, zP, XQ, YQ, XP, YP, Xlft, Xrght, Ylft, Yrght;

    for (P = 0; P < NVERTEX; ++P) {
        pvertex = vertex + P;
        ptr = pvertex->connect;
        if (ptr == NULL) continue;

        xP = pvertex->x; yP = pvertex->y; zP = pvertex->z;
        XP = xP / zP; YP = yP / zP;

        for (iconnect = 1; iconnect <= *ptr; ++iconnect) {
            Q = *(ptr + iconnect);
            pvertex = vertex + Q;
            xQ = pvertex->x; yQ = pvertex->y; zQ = pvertex->z;

//          drawsegment(xP, yP, zP, xQ, yQ, zQ); //!

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
            if (fabs(denom) <= eps)
                denom = eps;

            slope = (Yrght - Ylft) / denom;
            jbot = jtop = ywhole(Ylft);

            for (ipix = ipixleft; ipix <= ipixright; ++ipix) {
                if (ipix == ipixright)
                    jl = ywhole(Yrght);
                else
                    jl = ywhole(Ylft + (xreal(ipix + 1) - Xlft) * slope);

                lower[ipix] = min2(jbot,jl); jbot = jl;
                upper[ipix] = max2(jtop,jl); jtop = jl;
            }

            ntrset = 0;
            for (ipix = ipixleft; ipix <= ipixright; ++ipix)
                for (jpix = lower[ipix]; jpix <= upper[ipix]; ++jpix) {
                    pointer = &(ScreeN[ipix][jpix]);
                    pnode = pointer->start;
                    while (pnode != NULL) {
                        trnr = pnode->jtr;
                        /* tnrn will be stored only if it is not yet present in array trset(triangle set)*/
                        trset[ntrset] = trnr; /* sentinel*/

                        jtr = 0;
                        while (trset[jtr] != trnr)
                            jtr++
;
                        if (jtr == ntrset) {
                            ntrset++;
                            if (ntrset == NNTRSET) error("Triangle set overflow\n");
                        }
                        pnode = pnode->next;
                    }
                }
           drawsegment(xP, yP, zP, xQ, yQ, zQ);
           //linesegment(xP, yP, zP, xQ, yQ, zQ, 0);
        }
    }
}

void read_and_draw(double rho, double theta, double phi, double Xvp_min, double Xvp_max, double Yvp_min, double Yvp_max)
{
    /* Иниц. экранную матрицу, матрицу поворота, массив вершин, область вывода */
    extra_init(rho, theta, phi, Xvp_min, Xvp_max, Yvp_min, Yvp_max);

    /* Читаем вершины, преобразуем и записываем в массив vertex */
    read_vertex();

    /* Вычисляем значения для масштабирования */
    scale_calc(Xvp_min, Xvp_max, Yvp_min, Yvp_max);

    /* Читаем из файла грани обьекта и разбиваем на полигоны*/
    read_faces();

    /* Добавляем ближайшие треугольники в экранный список*/
    tr_add();

    /* Проверяем видимость прямых и рисуем их */
    line_and_pix();

}



void skipbl()
{
    char ch;
    do
        ch = getc(fpin);
    while
        (isspace(ch));
    ungetc(ch, fpin);
}


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

/*void add_linesegment(int P, int Q)
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
    vertex[P].connect = ptr = (int *)realloc(ptr, (n + 1) * isize);
    if (ptr == NULL) error("Memory alloc error \n");
    *(ptr + n) = Q; *ptr = n;
}

int counter_clock(int i0, int i1, int i2, double *pdist, int code)
{
    /*
        code = 0 : compute orientation
        code = 1 : compute a,b,c,h, store the first triangle
        code > 1 : check if next triangle is coplanar; store it
    /

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

    if (ntr = NTRIANGLE) error("Too many triangles");
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
	if (ipixleft != ipixright) slope = (Yright[i] - Yleft[i]) / denom;
	j_old = ywhole(Yleft[i]);
	for (ipix = ipixleft; ipix <= ipixright; ++ipix) {
	    if (ipix = ipixright)
                jl = xwhole(Yright[i]);
            else jl = xwhole(Yleft[i] + (xreal(ipix + 1) - Xleft[i]) * slope);
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
            pointer = &(Screen[ipix][jpix]);
            if (jpix > low[ipix] && jpix < up[ipix]) {
                xR = Xmin + (ipix + 0.5) * deltaX;
                yR = Ymin + (jpix + 0.5) * deltaY;
                denom = a * xR + b * yR + c * d;
                dist = fabs(denom) > eps ? h * sqrt(xR*xR + yR*yR + 1.) / denom : BIG;

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

*/

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
     //while (k < ntr) {
	ptriangle = triangle + k;
        a = ptriangle->a; b = ptriangle->b;  c = ptriangle->c;
        h = ptriangle->h;

       /* TEST #1*/

        hP = a * xP + b * yP + c * zP;
        hQ = a * xQ + b * yQ + c * zQ;
	eps1 = eps + eps * h;
        if (hP - h <= eps && hQ - h <= eps) {
            ++k;          /* PQ is not behind ABC*/
            continue;
        }

        /* TEST #2*/

        K1 = yP * zQ - yQ * zP; K2 = zP * xQ - zQ * xP; K3 = xP * yQ - xQ * yP;
        A = ptriangle->A;  B = ptriangle->B; C = ptriangle->C;
	pvertex = vertex + A;
        xA = pvertex->x; yA = pvertex->y; zA = pvertex->z;
	pvertex = vertex + B;
        xB = pvertex->x; yB = pvertex->y; zB = pvertex->z;
	pvertex = vertex + C;
        xC = pvertex->x; yC = pvertex->y; zC = pvertex->z;

        dA = K1 * xA + K2 * yA + K3 * zA;
        dB = K1 * xB + K2 * yB + K3 * zB;
        dC = K1 * xC + K2 * yC + K3 * zC;

        eA = dA > eps ? 1 : dA < meps ? -1 : 0;
        eB = dB > eps ? 1 : dB < meps ? -1 : 0;
        eC = dC > eps ? 1 : dC < meps ? -1 : 0;

        sum = eA + eB + eC;
        if (abs(sum) >= 2) {
            ++k;
            continue;
        }

	/* TEST #3*/

        Poutside = Qoutside = 0; labmin = 1; labmax = 0;
        for (i = 0; i < 3; ++i) {
            C1 = yA * zB - yB * zA; C2 = zA * xB - zB * xA; C3 = xA * yB - xB * yA;
            Cpos = C1 * xC + C2 * yC + C3 * zC;
            Ppos = C1 * xP + C2 * yP + C3 * zP;
            Qpos = C1 * xQ + C2 * yQ + C3 * zQ;

            denom1 = Qpos - Ppos;
            if (Cpos > eps) {
                Pbeyond = Ppos < meps; Qbeyond = Qpos < meps;
                outside = (Pbeyond && (Qpos <= eps)) || (Qbeyond && (Ppos <= eps));
            } else if (Cpos < meps) {
                       Pbeyond = Ppos > eps; Qbeyond = Qpos > eps;
                       outside = (Pbeyond && (Qpos >= meps)) || (Qbeyond && (Ppos >= meps));
                   } else outside = 1;

            if (outside)
                break;
            lab = fabs(denom1) <= eps ? 1.e7 : -Ppos / denom1;
                //lab указывает на точку пересечения PQ с EAB

            Poutside |= Pbeyond; Qoutside |= Qbeyond;
            denom2 = dB - dA;
            mu = fabs(denom2) <= eps ? 1.e7 : -dA / denom2;
                //mu указывает на точку пересечения AB с EPQ
            if (mu >= meps && mu <= oneplus && lab >= meps && lab <= oneplus) {
                if (lab < labmin)
                    labmin = lab;
                if (lab > labmax)
                    labmax = lab;
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
        } //PQ невидим


        /* TEST #5*/

        r1 = xQ - xP; r2 = yQ - yP; r3 = zQ - zP;

        xmin = xP + labmin * r1; ymin = yP + labmin * r2; zmin = zP + labmin * r3;
        if (a * xmin + b * ymin + c * zmin - h < -eps1) {
            ++k;
            continue;
        }

        xmax = xP + labmax * r2; ymax = yP + labmax * r2; zmax = zP + labmax * r3;
        if (a * xmax + b * ymax + c * zmax - h < -eps1) {
            ++k;
            continue;
        }

        /* TEST #6*/

        if (Poutside || hP < h - eps1)
            linesegment(xP, yP, zP, xmin, ymin, zmin, k + 1);
        if (Qoutside || hQ < h - eps1)
            linesegment(xQ, yQ, zQ, xmax, ymax, zmax, k + 1);
        worktodo = 0;
        break;
    }

    if (worktodo) {
        moveto((int)(d * xP / zP + c1), (int)(d * yP / zP + c2));
        lineto((int)(d * xQ / zQ + c1), (int)(d * yQ / zQ + c2));
    }
}

void init_viewport(double *pXMIN, double *pXMAX, double *pYMIN, double *pYMAX)
{
    double XMIN, XMAX, YMIN, YMAX;
    double len = 5;
    printf("Give viewport boundaries XMIN, XMAX, YMIN, YMAX\n");
    scanf("%lf %lf %lf %lf", pXMIN, pXMAX, pYMIN, pYMAX);

    XInitThreads();
    srand(time(NULL));

}

void opengraph()
{
    int gd = DETECT, gm;
    initgraph(&gd, &gm, NULL);
    setcolor(RED);
//    setfillstyle(1, WHITE);
}
