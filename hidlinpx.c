/*
*
*    Триангуляция только выпуклых многоугольников!
*    Память выделяется только под несколько интов
*    в fvertex элементы начинаются с №1 (а не с 0)
*    Границы вывода заданы по дефолту, чтобы избежать возможных проблем с памятью
*    ФУНКЦИЮ RAST_TR НЕ ТРОГАТЬ! Она работает сейчас как надо!
*
*    TODO
*       1. Освещение по фонгу(? при возможности
*       2. "Интерфейс"
*	    2.1 Перемещение центра координат
*	    2.2 Восстановление исходного масштаба
*	    2.3 Текст(?)
*       3. в Функцию ReadFace реализовать так, чтобы она читала каждую последовательность f/f/f пока не встретит '/0', и уже с ней работала
*
*/
#include<stdio.h>
#include<math.h>
#include<ctype.h>
#include<graphics.h>
#include "test_header.h"
#include<X11/Xlib.h>
#include<time.h>
#include<SDL/SDL.h>
#include<unistd.h>

int ntr = 0, nvr, lightNOW = 0;
double v11, v12, v13, v21, v22, v23, v32, v33, v43, d, c1, c2,
       Xrange, Yrange, Xmin, Xmax, Ymin, Ymax,
       deltaX, deltaY, xO = 0, yO = 0, zO = 0, onetime = 1,
       Xvp_range, Yvp_range;
const double Xvp_min = 10, Xvp_max = 450, Yvp_min = 10, Yvp_max = 450;

double eps = 1.e-5, meps = -1.e-5, oneminus = 1 - 1.e-5, oneplus = 1 + 1.e-5;

int *zbuffer;

struct color Button1C = {255, 0, 0};
struct color Button2C = {255, 0, 0};
struct color Black = {255, 255, 255};

struct ver *vertex, *pvertex, *fvertex;
struct tr *triangle, *ptriangle;
struct v3D lightdir = {0,0,1};
FILE *fpin, *fpout;

void drawsegment(double xP, double yP, double zP, double xQ, double yQ, double zQ);
void polydiv(int clockwise);
void update();
void check();
void update();
void border_draw(double XMIN, double XMAX, double YMIN, double YMAX);
void drawtr();
void extra_init();
void calc_vertex();
void scale_calc();
void calc_faces();
void tr_add();
void line_and_pix();
void read_and_draw(double rho, double theta, double phi, double scale);
void Rast_Tr(struct v2D A, struct v2D B, struct v2D C, struct color Cl);
void swap(double *x, double *y);
void swapS(struct v2D *A, struct v2D *B);
void normalize(struct ver vertex);
void ReadVertex(char *str);
void ReadFace(char *str);
void read_faces();
void calc_triangles();
void SPolyDiv(int *poly, int npoly);
void RememberVertex(int A, int B, int C, int Ntr);
void FinalizeMe();
void ReadNorm(char *str, int num);
void DrawButtons();
void RastRect(int x1, int y1, int x2, int y2, struct color Cl);


void DrawButtons()
{
    //Button 1
    DrawRect(Xvp_max, Yvp_min, Xvp_max + 20, Yvp_min + 20, Black, 2);
    RastRect(Xvp_max, Yvp_min, Xvp_max + 20, Yvp_min + 20, Button1C);

    //Button 2
    DrawRect(Xvp_max + 30, Yvp_min, Xvp_max + 50, Yvp_min + 20, Black, 2);
    RastRect(Xvp_max + 30, Yvp_min, Xvp_max + 50, Yvp_min + 20, Button2C);

    //screen
    DrawRect(Xvp_min, Yvp_min, Xvp_max, Yvp_max, Black, 1);
}

void swap(double *x, double *y)
{
    double aux;

    aux = *x;
    *x = *y;
    *y = aux;
}

void swapS(struct v2D *A, struct v2D *B)
{
    int xaux, yaux;
    xaux = A->x;
    A->x = B->x;
    B->x = xaux;

    yaux = A->y;
    A->y = B->y;
    B->y = yaux;

    double aux;
    aux = A->depth;
    A->depth = B->depth;
    B->depth = aux;
}

void RastRect(int x1, int y1, int x2, int y2, struct color Cl)
{
    while (x1 != x2 && y1 != y2) {
        x1 += 1; y1 += 1; x2 -= 1; y2 -= 1;
        DrawRect(x1, y1, x2, y2, Cl, 1);
    }
}

void Rast_Tr(struct v2D A, struct v2D B, struct v2D C, struct color Cl)
{
    if (A.y == B.y && A.y == C.y) return;

    if (A.y > B.y) swapS(&A,&B);
    if (A.y > C.y) swapS(&A,&C);
    if (B.y > C.y) swapS(&B,&C);

    int total_height = C.y - A.y, i, second_half, j, idx, segment_height;
    float alpha, beta, phi;
    double xI, yI, zI, xJ, yJ, zJ, Pz, Px, Py;

    for (i = 0; i <= total_height; ++i) {
        second_half = i > B.y - A.y || B.y == A.y;
        segment_height = second_half ? C.y - B.y : B.y - A.y;

        alpha = (float)(i) / total_height;
        beta = (float)(i - (second_half ? B.y - A.y : 0)) / segment_height;

        xI = alpha * C.x + A.x * (1 - alpha);
        yI = alpha * C.y + A.y * (1 - alpha);
	zI = alpha * C.depth + A.depth * (1 - alpha);

        xJ = second_half ? (beta * C.x + B.x * (1 - beta)) : beta * B.x + A.x * (1 - beta);
        yJ = second_half ? (beta * C.y + B.y * (1 - beta)) : beta * B.y + A.y * (1 - beta);
        zJ = second_half ? (beta * C.depth + B.depth * (1 - beta)) : beta * B.depth + A.depth * (1 - beta);

        if (xI > xJ) {
            swap(&xI, &xJ); swap(&yI, &yJ);
        }

        for (j = xI; j <= xJ; ++j) {
            if (j <= (Xvp_min) || j >= (Xvp_max) || (A.y + i) <= (Yvp_min) || (A.y + i) >= (Yvp_max) ) //будем считать, что окно вывода [50 350 50 350]
                continue;
	    phi = xJ == xI ? 1. : (float)(j-xI) / (float)(xJ-xI);
	    Px = xI + (xJ - xI) * phi;
	    Py = yI + (yJ - yI) * phi;
	    Pz = zI + (zJ - zI) * phi;
	    //idx = Px + 320 * Py;
	    idx = j + Yvp_range * (A.y + i);
	    if ((zbuffer + idx) == NULL) { printf("Buffer memory error\n"); continue; }
	    else if (*(zbuffer + idx) > Pz && *(zbuffer + idx) <= 20000) { //!*(zbuffer + idx) <= 20000 очень важный по-видимому момент, без него SegFault выдавало
                drawpix(j, A.y + i, Cl);
//		   printf("idx = %d, Pz = %lf, zbuffer + idx = %d \n", idx, Pz, *(zbuffer + idx));
		    *(zbuffer + idx) = Pz;
		//else {
		    //

		//}
	     //   drawpix(Px, Py, Cl);
	    }
        }
    }

//    update_render();
}

int main(int argc, char *argv[])
{
    int i, running = 1, autoROT = 0, rotate = 0, scaleON = 0,  xcountpos = 0, xcountneg = 0, ycountpos = 0, ycountneg = 0, yscalepos = 0, yscaleneg = 0;
    double rho = 3000, theta = 0, phi = 0, Ph = 0, T = 180, scale = 1, Xmaxe, Ymine, Xmine, Ymaxe;
    int dx, dy, x, y;

    init_viewport();
    init(1000,470,32);

    /* Открытие файла*/
    if (argc != 2 || (fpin = fopen(argv[1], "r")) == NULL)
        printf("Input file is not correcly specified");
    ntr = 0;
    nvr = 1;
    char string[100];
    int onemoretime = 1, nv = 0;
    while (!feof(fpin)) {
        fgets(string, sizeof(string), fpin);
        string[strlen(string) - 1] = '\0';
        switch (string[0]) {
            case 'v': if (string[1] == ' ') {
                          ReadVertex(string);
		      } else if (string[1] == 'n') {
			  ++nv; //счетчик вершин
		          ReadNorm(string, nv);
		      } else if (onemoretime && string[1] == 't') { //когда вершины для чтения кончатся, их нужно инициализировать для дальнейшего чтения
                          vertex = (struct ver*)malloc((nvr + 1) * sizeof(struct ver));
                          for (i = 1; i <= nvr; ++i) {
                              (vertex + i)->connect = (int *)malloc(100 * sizeof(int));
                              *((vertex + i)->connect) = 0;
                          }
		          onemoretime = 0;
		      }
                      break;
            case 'f': check(); ReadFace(string); break;
            default: continue; break;
        }
    }

    //read_vertex();
    //read_faces();

    read_and_draw(rho, theta, phi, scale);
    read_and_draw(rho, theta, phi, scale);

    int xscaleneg = 0, xscalepos = 0;
    SDL_Event event;
    while (running) {
        //Xmax = Xmax * d + c1; Xmin = Xmin * d + c1; Ymax = Ymax * d + c2; Ymin = Ymin * d + c2;
/*	if (lightNOW) {
	    T += 1; Ph += 1;
  //          update();
	    double factor = atan(1.0) / 45.0;
	    lightdir.x = sin(T * factor) * cos(Ph * factor);
    	    lightdir.y = sin(T * factor) * sin(Ph * factor);
            lightdir.z = cos(T * factor);

            read_and_draw(rho, theta, phi, scale);
	    sleep(0.9);
	    sleep(0.9);
	    sleep(0.9);
	}
*/	     while (autoROT) {
                theta += 1; phi += 1;
		read_and_draw(rho, theta, phi, scale);
                DrawButtons(Button1C, Button2C);
                if (SDL_PollEvent(&event)) {
                    if (event.type == SDL_MOUSEBUTTONDOWN && event.button.x > Xvp_max && event.button.x <= Xvp_max + 20 && event.button.y > Yvp_min && event.button.y < Yvp_min + 20) {
                        autoROT = 0;
                        Button1C.r = 255; Button1C.g = 0;
			read_and_draw(rho, theta, phi, scale);
			read_and_draw(rho, theta, phi, scale);

                    }
                }
	     }

        while (SDL_PollEvent(&event)) {
	    Xmaxe = Xmax * d + c1; Xmine = Xmin * d + c1; Ymaxe = Ymax * d + c2; Ymine = Ymin * d + c2;
	    switch (event.type) {
		case SDL_QUIT:
		    FinalizeMe();
                    running = 0;
		    break;
                case SDL_MOUSEBUTTONDOWN:
                    x = event.button.x;
                    y = event.button.y;

                    if (x > Xmine && x < Xmaxe && y > Ymine && y < Ymaxe  && event.button.button == SDL_BUTTON_LEFT)
                        rotate = 1;
		    if (x > Xvp_max + 30 && x < Xvp_max + 50 && y > Yvp_min && y < Yvp_min + 20) {
		//	lightNOW = !lightNOW;
			xO += 50; //yO += 5; zO += 5;
			Button2C.r = Button2C.r == 255 ? 0 : 255;
			Button2C.g = Button2C.g == 255 ? 0 : 255;
                        read_and_draw(rho, theta, phi, scale);
		    }

		    if (event.button.button == SDL_BUTTON_RIGHT && x < Xmaxe && x > Xmine && y > Ymine && y < Ymaxe) {
			scaleON = 1;
		    }
                    if (x > Xvp_max && x <= Xvp_max + 20 && y > Yvp_min && y < Yvp_min + 20) {
                        autoROT = 1;
                        Button1C.r = 0; Button1C.g = 255;
                    }

		    read_and_draw(rho, theta, phi, scale);

                    break;
                case SDL_MOUSEBUTTONUP:
                    rotate = 0;
		    scaleON = 0;
                    break;
                case SDL_MOUSEMOTION:
                    if (scaleON && event.motion.x < Xmaxe && event.motion.y < Ymaxe) {
                        dy = event.motion.yrel;
			dx = event.motion.xrel;

                        if (dy < 0) ++yscaleneg; if (dy > 0) ++yscalepos;
			if (dx < 0) ++xscaleneg; if (dx > 0) ++xscalepos;

                        if (yscalepos == 2) {
                            rho += 30; read_and_draw(rho, theta, phi, scale);
			    yscalepos = 0;
                        }

                        if (yscaleneg == 2) {
                            rho -= 30; read_and_draw(rho, theta, phi, scale);
			    yscaleneg = 0;
                        }
                    }
                    if (rotate && event.motion.x < Xmaxe && event.motion.y < Ymaxe) {
                        dx = event.motion.xrel; dy = event.motion.yrel;
                        if (dy > 0) ycountpos++; if (dy < 0) ycountneg++;
                        if (dx > 0) xcountpos++; if (dx < 0) xcountneg++;

                        if (ycountpos == 2) {
                            phi += 1; read_and_draw(rho, theta, phi, scale);
                            ycountpos = 0;
                        }
                        if (ycountneg == 2) {
                            phi -= 1; read_and_draw(rho, theta, phi, scale);
                            ycountneg = 0;
                        }
                        if (xcountpos == 2) {
                            theta +=1; read_and_draw(rho, theta, phi, scale);
                            xcountpos = 0;
                        }
                        if (xcountneg == 2) {
                            theta -= 1; read_and_draw(rho, theta, phi, scale);
                            xcountneg = 0;
                        }

                        break;
                    }
             }

        }
    }
}


void check()
{
   printf("checkpoint\n");
}

void update()
{
    int i;
    clearscreen();
    //if (zbuffer != NULL)
    //    free(zbuffer);
    for (i = 0; i < (Xvp_range + 1) * (Yvp_range + 1); ++i)
        if ( (zbuffer + i) != NULL)
	    *(zbuffer + i) = 20000;

}


void swapI(int *x, int *y)
{
    int aux;
    aux = *x;
    *x = *y;
    *y = aux;
}


void Rast_Guro(int n, struct v2D t0, struct v2D t1, struct v2D t2)
{
    struct color Cl, ClJ, ClI;
    double intensity, I1, I2, I3, IB, IA, IP;
    float gamma, Pz;
    int A = triangle[n].A, B = triangle[n].B, C = triangle[n].C, idx;

    if (t0.y > t1.y) {swapS(&t0,&t1); swapI(&A,&B);}
    if (t0.y > t2.y) {swapS(&t0,&t2); swapI(&A,&C);}
    if (t1.y > t2.y) {swapS(&t1,&t2); swapI(&B,&C);}

    int total_height = t2.y - t0.y;

    I1 = intensity = (vertex + A)->a * lightdir.x + (vertex + A)->b * lightdir.y + (vertex + A)->c * lightdir.z;
    struct color Cl1 = {255 * intensity, 255 * intensity, 255 * intensity};

    I2 = intensity = (vertex + B)->a * lightdir.x + (vertex + B)->b * lightdir.y + (vertex + B)->c * lightdir.z;
    struct color Cl2 = {255 * intensity, 255 * intensity, 255 * intensity};

    I3 = intensity = (vertex + C)->a * lightdir.x + (vertex + C)->b * lightdir.y + (vertex + C)->c * lightdir.z;
    struct color Cl3 = {255 * intensity, 255 * intensity, 255 * intensity};

    int second_half, segment_height;
   /* for (int i = 0.; i <= total_height; ++i) {
        second_half = i > t1.y - t0.y || t1.y == t0.y;
	IA = I1 * ()
    }*/
   for (int i = 0; i <= total_height; ++i) {
        second_half = i > t1.y - t0.y || t1.y == t0.y;
        segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;

        float alpha = (float)(i) / total_height;
        float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;

        double xI = alpha * t2.x + t0.x * (1 - alpha);
        double yI = alpha * t2.y + t0.y * (1 - alpha);
	double zI = alpha * t2.depth + t0.depth * (1 - alpha);

        double xJ = second_half ? (beta * t2.x + t1.x * (1 - beta)) : beta * t1.x + t0.x * (1 - beta);
        double yJ = second_half ? (beta * t2.y + t1.y * (1 - beta)) : beta * t1.y + t0.y * (1 - beta);
	double zJ = second_half ? (beta * t2.depth + t1.depth * (1 - beta)) : beta * t1.depth + t0.depth * (1 - beta);

        if (xI > xJ) {
            swap(&xI, &xJ); swap(&yI, &yJ);
        }
        IB = second_half ? I1 * (yJ - t1.y) / (t2.y - t1.y) + I2 * (t2.y - yJ) / (t2.y - t1.y) : I1 * (yJ - t0.y) / (t1.y - t0.y) + I2 * (t1.y - yJ) / (t1.y - t0.y);
        IA = I1 * (yI - t0.y) / (t2.y - t0.y) + I2 * (t2.y - yI) / (t2.y - t0.y);

        for (int j = xI; j <= xJ; ++j) {
            if (j <= Xvp_min || j >= Xvp_max || t0.y + i <= Yvp_min || t0.y + i >= Yvp_max) //будем считать, что окно вывода [50 350 50 350]
                continue;
            IP = xJ == xI ? IA : IA * (float)(xJ - j) / (float)(xJ - xI) + IB * (float)(j - xI) / (float)(xJ - xI);
            Cl.r = 255 * IP; Cl.g = 255 * IP; Cl.b = 255 * IP;

            gamma = xJ == xI ? 1. : (float)(j - xI) / (float)(xJ - xI);
            Pz = zI + (zJ - zI) * gamma;
            idx = j + Yvp_range * (t0.y + i);

            if (*(zbuffer + idx) > Pz) {
                drawpix(j, t0.y + i, Cl);
                *(zbuffer + idx) = Pz;
            }

        }
       /* for (int j = xI; j <= xJ; ++j) {
            if (j <= Xvp_min || j >= Xvp_max || t0.y + i <= Yvp_min || t0.y + i >= Yvp_max) //будем считать, что окно вывода [50 350 50 350]
                continue;

            gamma = xJ == xI ? 1. : (float)(j - xI) / (float)(xJ - xI);
            Cl.r = ClI.r * (1 - gamma) + ClJ.r * gamma;
            Cl.g = ClI.g * (1 - gamma) + ClJ.g * gamma;
            Cl.b = ClI.b * (1 - gamma) + ClJ.b * gamma;
	    Pz = zI + (zJ - zI) * gamma;
	    idx = j + Yvp_range * (t0.y + i);

	    if (*(zbuffer + idx) > Pz && Cl.r > 0 && Cl.g > 0 && Cl.b > 0) {
	        drawpix(j, t0.y + i, Cl);
		*(zbuffer + idx) = Pz;
	    }
	    //if (Cl.r > 0 && Cl.g > 0 && Cl.b > 0)
            //    drawpix(j, t0.y + i, Cl);
        }*/
    }

}
void drawtr()
{
    int i, XA, YA, XB, YB, XC, YC, A, B, C;
    double a, b, c, h, intensity, depth;

    for (i = 0; i < ntr; ++i) {
	a = triangle[i].a; b = triangle[i].b; c = triangle[i].c; h = triangle[i].h;
        A = triangle[i].A; B = triangle[i].B; C = triangle[i].C;

        intensity = lightdir.x * a + lightdir.y * b + lightdir.z * c;
        if (h > 0 && intensity > 0) {
        struct color Cl = {255 * intensity, 255 * intensity, 255 * intensity};

        XA = (int)( (vertex + A)->x * d / (vertex + A)->z + c1);
	YA = (int)( (vertex + A)->y * d / (vertex + A)->z + c2);
	XB = (int)( (vertex + B)->x * d / (vertex + B)->z + c1);
	YB = (int)( (vertex + B)->y * d / (vertex + B)->z + c2);
	XC = (int)( (vertex + C)->x * d / (vertex + C)->z + c1);
	YC = (int)( (vertex + C)->y * d / (vertex + C)->z + c2);
        int outA = 0, outB = 0, outC = 0;

        if ((XA <= Xvp_min || XA >= Xvp_max) || (YA <= Yvp_min || YA >= Yvp_max)) outA = 1;
        if ((XB <= Xvp_min || XB >= Xvp_max) || (YB <= Yvp_min || YB >= Yvp_max)) outB = 1;
        if ((XC <= Xvp_min || XC >= Xvp_max) || (YC <= Yvp_min || YC >= Yvp_max)) outC = 1;
        if (outA && outB && outC) continue;

        depth = (vertex + A)->z; struct v2D A = {XA, YA, depth};

        depth = (vertex + B)->z; struct v2D B = {XB, YB, depth};

        depth = (vertex + C)->z; struct v2D C = {XC, YC, depth};

        if (!lightNOW) Rast_Tr(A,B,C,Cl);
//	check();
 //       printf("Triangle (%d,%d)(%d, %d)(%d , %d) has been drown!\n", XA, YA, XB, YB, XC, YC);
	if (lightNOW) Rast_Guro(i, A, B, C);
        }
    }

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

void clear(char *str)
{
    int i;
    for (i = 0; i < strlen(str); ++i)
        str[i] = '\0';
}

void ReadFace(char *str)
{
    char numA[10], numB[10], numC[10], numD[10], num[10], part[30][150];
    int j = 0, i = 2, k = 0, A, B, C, D, poly[30];
    int IfNoSpaceAtTheEndOfTheLine = 1;

    while (1 == 1) {
        if (str[i] == ' ') {
            part[k][j] = '\0';
            if (str[i + 1] == '\0') break;
            ++i; ++k; j = 0;
            continue;
        } else if (str[i] == '\0') {
            part[k][j] = '\0';
            break;
        }
        part[k][j] = str[i];
        ++j; ++i;
    }

    for (i = 0; i <= k; ++i) {
        j = 0;
        while (part[i][j] != '/') {
            num[j] = part[i][j];
            ++j;
        }
        poly[i] = atoi(num);
//        clear(num);
	while (j >= 0) {
	    num[j] = '\0';
	    --j;
        }
    }

    SPolyDiv(poly, k);
}


/*void read_faces()
{
    int i;
    char ch;

//	triangle = (struct tr *)malloc(NTRIANGLE * sizeof(struct tr));

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


	simple_polydiv();
    }

    fclose(fpin);
}*/

void SPolyDiv(int *poly, int npoly)
{
    int i0 = 0, i1 = 1, i2 = 2, A, B, C;

    while(i2 != npoly + 1) {
        A = poly[i0], B = poly[i1], C = poly[i2];

        RememberVertex(A, B, C, ntr);

        if (triangle == NULL)
            triangle = (struct tr*)malloc(sizeof(struct tr));
        else
            triangle = (struct tr*)realloc(triangle, (ntr + 1) * sizeof(struct tr));

        (triangle + ntr)->A = A; (triangle + ntr)->B = B; (triangle + ntr)->C = C;
        ++ntr; ++i1; ++i2;

        printf("Треугольник №%d / 2492 добавлен успешно!\n", ntr);
    }
}


void extra_init()
{

    int i;
    zbuffer = (int *)malloc((Xvp_range + 1)* (Yvp_range + 1)* sizeof(int));
    for (i = 0; i < (Xvp_range + 1) * (Yvp_range + 1); ++i)
            *(zbuffer + i) = 20000;
}

void ReadVertex(char *str)
{
//    char ch;
    char numx[15], numy[15], numz[15];

    int i = 2, j = 0;
    while (str[i] != ' ') {
        numx[j] = str[i];
        ++i; ++j;
    }
    double x = atof(numx) * 100;

    ++i; j = 0;
    while (str[i] != ' ') {
        numy[j] = str[i];
        ++i; ++j;
    }
    double y = atof(numy) * 100;

    ++i; j = 0;
    while (str[i] != '\0') {
        numz[j] = str[i];
       ++i; ++j;
    }
    double z = atof(numz) * 100;

/*    fscanf(fpin, "%lf %lf %lf", &xO, &yO, &zO);
    while(skipbl(), ch = getc(fpin), ch != 'F' && ch != 'f') {
        ungetc(ch,fpin);
        fscanf(fpin, "%d %lf %lf %lf", &i, &x, &y, &z);
        if (i < 0 || i >= NVERTEX )
            error("Неверный номер вершины \n");*/
//check();
    if (fvertex != NULL)
        fvertex = (struct ver *)realloc(fvertex, (nvr + 1) * sizeof(struct ver));
    else fvertex = (struct ver *)malloc(2 * sizeof(struct ver));
//check();
        (fvertex + nvr)->x = x; (fvertex + nvr)->y = y; (fvertex + nvr)->z = z;
        nvr++;
   // }

   /* for (i = 0; i < nvr; ++i) {
	vertex[i].connect = (int *)malloc(100 * sizeof(int));
	*(vertex[i].connect) = 0;
    }*/
}

void ReadNorm(char *str, int num)
{
    char numx[15], numy[15], numz[15];

    int i = 4, j = 0;
    while (str[i] != ' ') {
        numx[j] = str[i];
        ++i; ++j;
    }
    double x = atof(numx);

    ++i; j = 0;
    while (str[i] != ' ') {
        numy[j] = str[i];
        ++i; ++j;
    }
    double y = atof(numy);

    ++i; j = 0;
    while (str[i] != '\0') {
        numz[j] = str[i];
       ++i; ++j;
    }
    double z = atof(numz);

    (fvertex + num)->a = x; (fvertex + num)->b = y; (fvertex + num)->c = z;

}

void calc_vertex()
{
    int i;
    double x, y, z, xe, ye, ze, X, Y;

    Xmin = Ymin = BIG; Xmax = Ymax = -BIG;
    for (i = 1; i <= nvr; ++i) {
        x = (fvertex + i)->x; y = (fvertex + i)->y; z = (fvertex + i)->z;

        //преобразуем координаты (через матрицу поворота)
        viewing(x - xO, y - yO, z - zO, &xe, &ye, &ze);
//        if (ze <= eps)
//            continue;

        //находим max и min для масштабирования в наше окно вывода
        X = xe / ze; Y = ye / ze;
        if (X < Xmin) Xmin = X; if (X > Xmax) Xmax = X;
        if (Y < Ymin) Ymin = Y; if (Y > Ymax) Ymax = Y;

        (vertex + i)->x = xe; (vertex + i)->y = ye; (vertex + i)->z = ze;
//	printf("%lf %lf %lf\n", xe, ye, ze);
    }
 //   printf("endofinput\n");
}

void calc_ver2()
{
    int i, j, *ptr, Ntr;
    double a, b, c, r;
    for (i = 1; i <= nvr; ++i) {
	a = 0; b = 0; c = 0;
        ptr = (vertex + i)->connect;

        Ntr = *ptr;

        for (j = 1; j <= Ntr; ++j) {
            a += (triangle + *(ptr + j) )->a;
            b += (triangle + *(ptr + j) )->b;
            c += (triangle + *(ptr + j) )->c;
        }

        r = sqrt(a*a + b*b + c*c);
        a /= r; b /= r; c /= r;

        (vertex + i)->a = a; (vertex + i)->b = b; (vertex + i)->c = c;
    }

}

void calc_triangles()
{
    double xA, yA, zA, xB, yB, zB, xC, yC, zC, a, b, c, h, r;
    int i, A, B, C;
    for (i = 0; i < ntr; ++i) {
        A = (triangle + i)->A; B = (triangle + i)->B; C = (triangle + i)->C;

	xA = (vertex + A)->x; yA = (vertex + A)->y; zA = (vertex + A)->z;
        xB = (vertex + B)->x; yB = (vertex + B)->y; zB = (vertex + B)->z;
        xC = (vertex + C)->x; yC = (vertex + C)->y; zC = (vertex + C)->z;

        a = yA * (zB - zC) - yB * (zA - zC) + yC * (zA - zB);
        b = -(xA * (zB - zC) - xB * (zA - zC) + xC * (zA - zB));
        c =  xA * (yB - yC) - xB * (yA - yC) + xC * (yA - yB);
	h = xA * (yB * zC - yC * zB) - xB * (yA * zC - yC * zA) + xC * (yA * zB - yB * zA);

        r = sqrt(a*a + b*b + c*c); r = r == 0 ? eps : r;
        a /= r; b /= r; c /= r;
	h /= r;

        (triangle + i)->a = a; (triangle + i)->b = b; (triangle + i)->c = c;
	(triangle + i)->h = h;
    }
}

void scale_calc()
{
    double fx, fy, Xcentre, Ycentre, Xvp_centre, Yvp_centre;

    Xrange = Xmax - Xmin; Yrange = Ymax - Ymin;
    Xvp_range = Xvp_max - Xvp_min; Yvp_range = Yvp_max - Yvp_min;
    fx = Xvp_range / Xrange; fy = Yvp_range / Yrange;
    d = (fx < fy) ? fx : fy;
    Xcentre = (Xmin + Xmax) / 2; Ycentre = (Ymax + Ymin) / 2;
    Xvp_centre = (Xvp_min + Xvp_max) / 2; Yvp_centre = (Yvp_min + Yvp_max) / 2;
    c1 = Xvp_centre - d * Xcentre; c2 = Yvp_centre - d * Ycentre;
}

void read_and_draw(double rho, double theta, double phi, double scale)
{

    update();
//check();
    coeff(rho, theta, phi, scale);

    calc_vertex();
//check();
    if (onetime) {
        scale_calc();
	extra_init();
    }
    onetime = 0;

//check();
     //coeff(rho, theta, phi, scale);
    /* Иниц. z-буфер */
   // 
//check();
    /* Читаем вершины, преобразуем и записываем в массив vertex */
//    calc_vertex();
//check();
    /* Вычисляем значения для масштабирования (один раз)*/

//check();
    calc_triangles();
//check();
    /* Считаем нормали к вершинам*/
    calc_ver2();
//check();
//    struct color Cl = {174, 32, 14};
//    rectangle(Xmaxe, Ymine, Xmaxe + 50, Ymine + 50, Cl);
    drawtr();
    DrawButtons();

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


void RememberVertex(int A, int B, int C, int Ntr)
{
   int *ptr, n;

   ptr = (vertex + A)->connect; n = *ptr; ++n;
//   ptr = (int *)realloc(ptr, (n + 1) * sizeof(int));
   *(ptr + n) = Ntr; *ptr = n;
//check();
   ptr = (vertex + B)->connect; n = *ptr; ++n;
  // ptr = (int *)realloc(ptr, (n + 1) * sizeof(int));
   *(ptr + n) = Ntr; *ptr = n;
//check();
   ptr = (vertex + C)->connect; n = *ptr; ++n;
  // ptr = (int *)realloc(ptr, (n + 1) * sizeof(int));
   *(ptr + n) = Ntr; *ptr = n;
//check();
}

void error(char *string)
{
    closegraph();
    printf("%s \n", string);
}

void coeff(double rho, double theta, double phi, double scale)
{
    double th, ph, costh, sinth, cosph, sinph, factor;

    factor = atan(1.0) / 45.0;
    th = theta * factor; ph = phi * factor;
    costh = cos(th); sinth = sin(th);
    cosph = cos(ph); sinph = sin(ph);

    v11 = -sinth * scale; v12 = -cosph * costh * scale; v13 = -sinph * costh * scale;
    v21 = costh * scale;  v22 = -cosph * sinth * scale; v23 = -sinph * sinth * scale;
                          v32 = sinph * scale;          v33 = -cosph * scale;
                                                        v43 = rho;

}

void viewing(double x, double y, double z, double *pxe, double *pye, double *pze)
{
    *pxe = (v11 * x + v21 * y);
    *pye = (v12 * x + v22 * y + v32 * z);
    *pze = (v13 * x + v23 * y + v33 * z + v43);

}


void init_viewport()
{
    //printf("Give viewport boundaries XMIN, XMAX, YMIN, YMAX\n");
    //scanf("%lf %lf %lf %lf", &Xvp_min, &Xvp_max, &Yvp_min, &Yvp_max);

    //XInitThreads();
    //srand(time(NULL));
}

void FinalizeMe()
{
    int i;

    free(zbuffer);
    free(triangle);
    for (i = 1; i <= nvr; ++i)
        free(vertex[i].connect);
    free(vertex);
    free(fvertex);
}
