/*
*
*    Триангуляция только выпуклых многоугольников!
*    Память выделяется только под несколько интов
*    в fvertex элементы начинаются с №1 (а не с 0)
*    Границы вывода заданы по дефолту, чтобы избежать возможных проблем с памятью
*    ФУНКЦИЮ RAST_TR НЕ ТРОГАТЬ! Она работает сейчас как надо!
*    Если в функции MoveObj() вершины начинать считать с 0, тоесть xO yO zO, то выходит(?) что у каждого обьекта своя система координат со своим центром
*		ну а если нет, то своя система только у главного обьекта, пиздец
*    Разобраться с геометрией
*
*    TODO
*       1. Освещение по гуро(? при возможности)
*       2. "Интерфейс"
*	    2.1 Перемещение центра координат
*	    2.2 Восстановление исходного масштаба
*	    2.3 Текст(?)
*       3. Чтение .obj файла в котором несколько обьектов (массив структур обьектов)
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

int lightNOW = 0;
double v11, v12, v13, v21, v22, v23, v32, v33, v43, d, c1, c2, xO = 0,// yO = 0, zO = 0
       onetime = 1, Xvp_range, Yvp_range;
const double Xvp_min = 10, Xvp_max = 450, Yvp_min = 10, Yvp_max = 450;
double eps = 1.e-5, meps = -1.e-5, oneminus = 1 - 1.e-5, oneplus = 1 + 1.e-5;

int *zbuffer;

struct color Button1C = {255, 0, 0};
struct color Button2C = {255, 0, 0};
struct color Black = {255, 255, 255};

//struct ver *vertex, *pvertex, *fvertex;
//struct tr *triangle, *ptriangle;
struct v3D lightdir = {0,0,1};

void update();
void check();
void DrawTr(struct tr* triangle, int ntr, struct ver* vertex);
void extra_init();
void CalcVertex(struct ver* fvertex, struct ver** vertex, int nvr, double *Range);
void extraCalcVertex(struct ver** vertex, int nvr, struct tr* triangle);
void scale_calc(double Xmin, double Xmax, double Ymin, double Ymax);
void CalcTr(struct ver* vertex, struct tr** triangle, int ntr);
void read_and_draw(struct ver** fvertex, struct ver** vertex, int nvr, struct tr** triangle, int ntr, double rho, double theta, double phi);
void Rast_Tr(struct v2D A, struct v2D B, struct v2D C, struct color Cl);
void swap(double *x, double *y);
void swapS(struct v2D *A, struct v2D *B);
void normalizeMe(struct v3D *vec);
void ReadVertex(char *str, struct ver **fvertex, int *nvr);
void ReadFace(char *str, struct tr** triangle, int *ntr);
void SPolyDiv(int *poly, int npoly, struct tr* triangle, int *ntr);
void RememberVertex(struct tr* triangle, int ntr, struct ver** vertex);
void FinalizeMe(struct ver** vertex, int nvr, struct ver** fvertex, struct tr** triangle);
//void ReadNorm(char *str, int num);
void DrawButtons();
void RastRect(int x1, int y1, int x2, int y2, struct color Cl);
void VerInit(struct ver** vertex, int nvr);
void Rast_Texture(char* filename, struct v2D A, struct v2D B, struct v2D C);
int Proection(char flag, struct ver* vertex, int num);
int isOutside(struct v2D A, struct v2D B, struct v2D C);

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

void Rast_Texture(char* filename, struct v2D A, struct v2D B, struct v2D C)
{
    SDL_Surface *image;
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
            idx = j + Yvp_range * (A.y + i);
            if ((zbuffer + idx) == NULL) { printf("Buffer memory error\n"); continue; }
            else if (*(zbuffer + idx) > Pz && *(zbuffer + idx) <= 20000) { //!*(zbuffer + idx) <= 20000 очень важный по-видимому момент, без него SegFault выдает
                LoadTexture(image, j, A.y + i);
                *(zbuffer + idx) = Pz;
            }
        }
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
		*(zbuffer + idx) = Pz;
	    }
        }
    }

}

void VerInit(struct ver **vertex, int nvr)
{   int i;
    *vertex = (struct ver*)malloc((nvr + 1) * sizeof(struct ver));
    for (i = 1; i <= nvr; ++i) {
        (*vertex + i)->connect = (int *)malloc(100 * sizeof(int));
        *((*vertex + i)->connect) = 0;
    }

}

void moveObj(struct ver* vertex, int nvr, int dx, int dy, int dz)
{
    int i;
    for (i = 1; i <= nvr; ++i) {
        vertex[i].x += dx;
	vertex[i].y += dy;
	vertex[i].z += dz;
    }
}

void rotateObj(struct ver* fvertex, int nvr, double theta, double phi)
{
    int i;
    double xe, ye, ze;
    coeff(1, theta, phi); //rho = 3000
    for (i = 1; i <= nvr; ++i) {
        viewing(fvertex[i].x, 0, fvertex[i].y, 0, fvertex[i].z, 0, &xe, &ye, &ze);
//	printf("Исходные координаты: %lf %lf %lf \n", fvertex[i].x, fvertex[i].y, fvertex[i].z);
// 	printf("Преобразованные: %lf %lf %lf \n", xe, ye, ze);
	fvertex[i].x = xe; fvertex[i].y = ye; fvertex[i].z = ze;
    }
}

void ReadFromObj(char *filename, struct ver **fvertex, struct ver** vertex, int *nvr, struct tr** triangle, int* ntr)
{
    FILE *fpin;
    char string[100];
    int onemoretime = 1, nv = 0;

    fpin = fopen(filename, "r");
    while (!feof(fpin)) {
        fgets(string, sizeof(string), fpin);
        string[strlen(string) - 1] = '\0';
        switch (string[0]) {
            case 'v': if (string[1] == ' ') {
                          ReadVertex(string, fvertex, nvr);
                      } else if (onemoretime) { //когда вершины для чтения кончатся, их нужно инициализировать для дальнейшего чтения
                          VerInit(vertex, *nvr);
                          onemoretime = 0;
                      }
                      break;
            case 'f': ReadFace(string, triangle, ntr); break;
            default: continue; break;
        }
    }
    fclose(fpin);
    RememberVertex(*triangle, *ntr, vertex);

}

void GameOn(struct ver **fvertex, struct ver **vertex, int nvr, struct tr** triangle, int ntr)
{
    int i, j;
    float rho = 3000, theta, phi;
    printf("1. %lf %lf %lf \n", (*fvertex)[1].x, (*fvertex)[1].y, (*fvertex)[1].z);
   //Загрузка обьектов //
    struct ver* v_plane = NULL, *fv_plane = NULL; struct tr* plane = NULL; int PlaneNtr = 0, PlaneNvr = 1;
    ReadFromObj("plane.obj", &fv_plane, &v_plane, &PlaneNvr, &plane, &PlaneNtr);

    struct tr *cube[3]; struct ver *v_cube[3], *fv_cube[3]; int CubeNtr[3] = {0,0,0}, CubeNvr[3] = {1, 1, 1};

    for (i = 0; i < 3; ++i)
        ReadFromObj("cube.obj", &fv_cube[i], &v_cube[i], &CubeNvr[i], &cube[i], &CubeNtr[i]);

    for (i = 0; i < 3; ++i ) {
        moveObj(fv_cube[i], CubeNvr[i], 0, i * 10, 200);
    }

/*    for (phi = 0, theta = 0; ;theta += 1, phi += 1, d = d >= 400 ? d - 25 : 399) {
        update();
	moveObj(*fvertex, nvr, dx, 0, 0);
	    counter += dx;
	    if (counter >= 2300) dx = 0;
        read_and_draw(fvertex, vertex, nvr, triangle, ntr, rho, theta, phi);
        read_and_draw(&fv_plane, &v_plane, PlaneNvr, &plane, PlaneNtr, rho, theta, phi);
        for (j = 0; j < 3; ++j)
            read_and_draw(&fv_cube[j], &v_cube[j], CubeNvr[j], &cube[j], CubeNtr[j], rho, theta, phi);
	if (dx == 0 && (int)theta % 90 == 0 && (int)phi % 360 == 0)
	    break;
    }
    getchar();*/
    //moveObj(*fvertex, nvr, 2300, -100, 0);
    d = 300;
    theta = 0; phi = 100;
    printf("1. %lf %lf %lf \n", (*fvertex)[1].x, (*fvertex)[1].y, (*fvertex)[1].z);
    rotateObj(*fvertex, nvr, 90, 90);
    rotateObj(*fvertex, nvr, 180, 180);
    printf("2. %lf %lf %lf \n", (*fvertex)[1].x, (*fvertex)[1].y, (*fvertex)[1].z);
    for (int i = 0; i < 60; i++) {
	sleep(0.9);
        update();
        moveObj(*fvertex, nvr, 30, 0 ,0);
	//printf("#. %lf %lf %lf \n", (*fvertex)[1].x, (*fvertex)[1].y, (*fvertex)[1].z);
        read_and_draw(fvertex, vertex, nvr, triangle, ntr, rho, theta, phi);
        read_and_draw(&fv_plane, &v_plane, PlaneNvr, &plane, PlaneNtr, rho, theta, phi);
        //for (j = 0; j < 3; ++j)
        //    read_and_draw(&fv_cube[j], &v_cube[j], CubeNvr[j], &cube[j], CubeNtr[j], rho, theta, phi);
    }

    int run = 1, py = 0, keyHeld[323]= {0}; double phi_o = 0;
    SDL_Event event;
    printf("%lf %lf %lf \n", fv_cube[0][2].x, fv_cube[0][2].y, fv_cube[0][2].z);
    while (run) {
        update();
        for (i = 0; i < 3; ++i) {
            //moveObj(fv_cube[i], CubeNvr[i], 0, 10 * (i + 1) * 0.22, - 10 * (i + 1));
//	    rotateObj(*fvertex, nvr, 0, );
	    moveObj(fv_cube[i], CubeNvr[i], 10 * (i + 1), 0, 0);
            read_and_draw(&fv_plane, &v_plane, PlaneNvr, &plane, PlaneNtr, rho, theta, phi); //plane
            read_and_draw(&fv_cube[i], &v_cube[i], CubeNvr[i], &cube[i], CubeNtr[i], rho, theta, phi); //enemies
	    read_and_draw(fvertex, vertex, nvr, triangle, ntr, rho, theta, phi);//main character
	    if (fv_cube[i][1].x > 2300) {
		printf("%d -ый куб : %lf %lf %lf \n", i, fv_cube[i][2].x, fv_cube[i][2].y, fv_cube[i][2].z);
		printf("Координаты обьекта: %lf %lf %lf\n", (*fvertex)[2].x, (*fvertex)[2].y, (*fvertex)[2].z);
	        py += 120;
	        py = (fv_cube[i][2].y + py) > 500 ? py = -(fv_cube[i][2].y + 500) : py;
	        moveObj(fv_cube[i], CubeNvr[i], -4600, py, 0);
	    }
        }
        if (SDL_PollEvent(&event)) {
            if (event.type == SDL_MOUSEBUTTONDOWN && event.button.x > Xvp_max)
	        run = 0;
            if (event.type == SDL_KEYDOWN) {
                keyHeld[event.key.keysym.sym] = 1;
            }
	    if (event.type == SDL_KEYUP) {
                keyHeld[event.key.keysym.sym] = 0;
	    }
        }

        if (keyHeld[SDLK_LEFT])
            moveObj(*fvertex, nvr, 0, -10, 0);
        if (keyHeld[SDLK_RIGHT])
	    moveObj(*fvertex, nvr, 0, 10, 0);
    }

    return;

}
int main(int argc, char *argv[])
{
    int i, running = 1, autoROT = 0, rotate = 0, scaleON = 0,  xcountpos = 0, xcountneg = 0, ycountpos = 0, ycountneg = 0, yscalepos = 0, yscaleneg = 0;
    double rho = 3000, theta = 0, phi = 0, Ph = 0, T = 180, scale = 1, Xmaxe, Ymine, Xmine, Ymaxe;
    int dx, dy, x, y;

    int keyHeld[323]= {0};

    struct tr *test_triangle;
    struct ver *test_vertex, *test_fvertex;
    int test_nvr = 1, test_ntr = 0;

    init(1000,470,32);

    /* Открытие файла*/
    ReadFromObj(argv[1], &test_fvertex, &test_vertex, &test_nvr, &test_triangle, &test_ntr);
    /* Чтение кубов */
   // struct tr *cube[4];
   // struct ver *vercube[4], *fvercube[4];
  ///  int ntr[4] = {0,0,0, 0}, nvr[4] = {1, 1, 1, 1};

   // for (i = 0; i < 4; ++i)
    //    ReadFromObj("cube.obj", &fvercube[i], &vercube[i], &nvr[i], &cube[i], &ntr[i]);

   // for (i = 0; i < 3; ++i)
//	moveObj(fvercube[i], nvr[i], 0, i * 10, 0);
 //   moveObj(fvercube[3], nvr[3], 4600, 100, 0);

    read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);
//    read_and_draw(&fvercube[0], &vercube[0], nvr[0], &cube[0], ntr[0], rho, theta, phi);


    int xscaleneg = 0, xscalepos = 0, ObjectMoves = 0, py = 0;
    SDL_Event event;
    while (running) {
      /*  while (ObjectMoves) {
	    update();
	    for (i = 0; i < 3; ++i) {
	        moveObj(fvercube[i], nvr[i], 10 * (i + 1), 0, 0);
	        read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);
                read_and_draw(&fvercube[i], &vercube[i], nvr[i], &cube[i], ntr[i], rho, theta, phi); //enemies
	        read_and_draw(&fvercube[3], &vercube[3], nvr[3], &cube[3], ntr[3], rho, theta, phi); //main character
	        if (fvercube[i][1].x >= 2300) {
	            py += 120;
	 	    py = (fvercube[i][2].y + py) > 500 ? py = -(fvercube[i][2].y + 500) : py;
		    moveObj(fvercube[i], nvr[i], -4600, py, 0);
	        }

	    }
	    if (SDL_PollEvent(&event)) {
		if (event.type == SDL_MOUSEBUTTONDOWN && event.button.x > Xvp_max)
		    ObjectMoves = 0;

		if (event.type == SDL_KEYDOWN) {
		    keyHeld[event.key.keysym.sym] = 1;
		}
		if (event.type == SDL_KEYUP) {
		    keyHeld[event.key.keysym.sym] = 0;
		}
	    }

            if (keyHeld[SDLK_LEFT])
		moveObj(fvercube[3], nvr[3], 0, -10, 0);
	    if (keyHeld[SDLK_RIGHT])
		moveObj(fvercube[3], nvr[3], 0, 10, 0);
	}
*/
        while (autoROT) {
            theta += 1; phi += 1;
	    update();
	    read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);
	    //read_and_draw(&fvercube[0], &vercube[0], nvr[0], &cube[0], ntr[0], rho, theta, phi);
            DrawButtons(Button1C, Button2C);
            if (SDL_PollEvent(&event)) {
                if (event.type == SDL_MOUSEBUTTONDOWN && event.button.x > Xvp_max && event.button.x <= Xvp_max + 20 && event.button.y > Yvp_min && event.button.y < Yvp_min + 20) {
                    autoROT = 0;
                    Button1C.r = 255; Button1C.g = 0;
	            read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);
	//	    read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);

                }
            }
        }

        while (SDL_PollEvent(&event)) {
//	    Xmaxe = Xmax * d + c1; Xmine = Xmin * d + c1; Ymaxe = Ymax * d + c2; Ymine = Ymin * d + c2;
	    switch (event.type) {
		case SDL_QUIT:
		    FinalizeMe(&test_vertex, test_nvr, &test_fvertex, &test_triangle);
		   // for (i = 0; i < 4; ++i)
		//	FinalizeMe(&vercube[i], nvr[i], &fvercube[i], &cube[i]);
                    running = 0;
		    break;
                case SDL_MOUSEBUTTONDOWN:
                    x = event.button.x;
                    y = event.button.y;

                    if (x > Xvp_min && x < Xvp_max && y > Yvp_min && y < Yvp_max  && event.button.button == SDL_BUTTON_LEFT)
                        rotate = 1;
		    if (x > Xvp_max + 30 && x < Xvp_max + 50 && y > Yvp_min && y < Yvp_min + 20) {
		//	lightNOW = !lightNOW;
			//xO += 50; //yO += 5; zO += 5;
			Button2C.r = Button2C.r == 255 ? 0 : 255;
			Button2C.g = Button2C.g == 255 ? 0 : 255;
			//ObjectMoves = 1;
			GameOn(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr);
                        read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);
		    }

		    if (event.button.button == SDL_BUTTON_RIGHT && x < Xvp_max && x > Xvp_min && y > Yvp_min && y < Yvp_max) {
			scaleON = 1;
		    }
                    if (x > Xvp_max && x <= Xvp_max + 20 && y > Yvp_min && y < Yvp_min + 20) {
                        autoROT = 1;
                        Button1C.r = 0; Button1C.g = 255;
                    }

		    read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);

                    break;
                case SDL_MOUSEBUTTONUP:
                    rotate = 0;
		    scaleON = 0;
                    break;
                case SDL_MOUSEMOTION:
                    if (scaleON && event.motion.x < Xvp_max && event.motion.y < Yvp_max) {
                        dy = event.motion.yrel;
			dx = event.motion.xrel;

                        if (dy < 0) ++yscaleneg; if (dy > 0) ++yscalepos;
			if (dx < 0) ++xscaleneg; if (dx > 0) ++xscalepos;

                        if (yscalepos == 2) {
                            d += 30;
			    yscalepos = 0;
                        }

                        if (yscaleneg == 2) {
                            d -= 30;
			    yscaleneg = 0;
                        }
			update();
			read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);
                    }
                    if (rotate && event.motion.x < Xvp_max && event.motion.y < Yvp_max) {
                        dx = event.motion.xrel;
			dy = event.motion.yrel;
                        if (dy > 0) ycountpos++; if (dy < 0) ycountneg++;
                        if (dx > 0) xcountpos++; if (dx < 0) xcountneg++;

                        if (ycountpos == 2) {
                            phi += 1;
                            ycountpos = 0;
                        }
                        if (ycountneg == 2) {
                            phi -= 1;
                            ycountneg = 0;
                        }
                        if (xcountpos == 2) {
                            theta +=1;
                            xcountpos = 0;
                        }
                        if (xcountneg == 2) {
                            theta -= 1;
                            xcountneg = 0;
                        }
		        update();
			printf("theta = %f phi = %f\n", theta, phi);
		        read_and_draw(&test_fvertex, &test_vertex, test_nvr, &test_triangle, test_ntr, rho, theta, phi);

                        //break;
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


void Rast_Guro(struct tr* triangle, struct ver* vertex, int n, struct v2D t0, struct v2D t1, struct v2D t2)
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

int Proection(char flag, struct ver* vertex, int num)
{
    switch (flag) {
	case 'x' : return (int)( vertex[num].x * d / vertex[num].z + c1); break;
        case 'y' : return (int)( vertex[num].y * d / vertex[num].z + c2); break;
	default : return 1;
    }
}

void normalizeMe(struct v3D *vec)
{
    double a, b, c, r;
    r = sqrt(sqr(vec->x) + sqr(vec->y) + sqr(vec->z));
    vec->x /= r; vec->y /= r; vec->z /= r;
}

int isOutside(struct v2D A, struct v2D B, struct v2D C)
{
    int outA = 0, outB = 0, outC = 0;

    if ((A.x <= Xvp_min || A.x >= Xvp_max) || (A.y <= Yvp_min || A.y >= Yvp_max)) outA = 1;
    if ((B.x <= Xvp_min || B.x >= Xvp_max) || (B.y <= Yvp_min || B.y >= Yvp_max)) outB = 1;
    if ((C.x <= Xvp_min || C.x >= Xvp_max) || (C.y <= Yvp_min || C.y >= Yvp_max)) outC = 1;

    if (outA && outB && outC) return 1;
    else return 0; 

}

void DrawTr(struct tr* triangle, int ntr, struct ver* vertex)
{
    int i, XA, YA, XB, YB, XC, YC, A, B, C, clA, clB, clC;
    double a, b, c, h, intensity, depth;
    normalizeMe(&lightdir);
    for (i = 0; i < ntr; ++i) {
	a = triangle[i].a; b = triangle[i].b; c = triangle[i].c; h = triangle[i].h;
        A = triangle[i].A; B = triangle[i].B; C = triangle[i].C;

        intensity = lightdir.x * a + lightdir.y * b + lightdir.z * c;
        if (h > 0 && intensity > 0) {

        //struct color Cl = {(int)(a * 1000) % 255 * intensity, (int)(b * 1000) % 255 * intensity, (int)(c * 1000) % 255 * intensity};
	//struct color Cl = {17, 16, 174};
        struct color Cl = {255 * intensity, 255 * intensity, 255 * intensity};
        XA = Proection('x', vertex, A); YA = Proection('y', vertex, A);
	XB = Proection('x', vertex, B); YB = Proection('y', vertex, B);
	XC = Proection('x', vertex, C); YC = Proection('y', vertex, C);

        struct v2D vA = {XA, YA, vertex[A].z};
        struct v2D vB = {XB, YB, vertex[B].z};
        struct v2D vC = {XC, YC, vertex[C].z};

        if (isOutside(vA, vB, vC)) continue;

        if (!lightNOW) Rast_Tr(vA,vB,vC,Cl);
//	Rast_Texture("grass2.bmp", A, B, C);
 //       printf("Triangle (%d,%d)(%d, %d)(%d , %d) has been drown!\n", XA, YA, XB, YB, XC, YC);
 //       Rast_Guro(triangle, vertex, i, A, B, C);
        }
    }

}

void clear(char *str)
{
    int i;
    for (i = 0; i < strlen(str); ++i)
        str[i] = '\0';
}

void ReadFace(char *str, struct tr **triangle, int *ntr)
{
    char numA[10], numB[10], numC[10], numD[10], num[10], part[30][150];
    int j = 0, i = 2, k = 0, poly[30];

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

	while (j >= 0) {
	    num[j] = '\0';
	    --j;
        }
    }

    //SPolyDiv(poly, k, triangle, ntr);
    int i0 = 0, i1 = 1, i2 = 2, A, B, C;

    while(i2 != k + 1) {
        A = poly[i0], B = poly[i1], C = poly[i2];

    //    RememberVertex(A, B, C, ntr);
        printf("%d\n", *ntr);
        if (*ntr == 0)
            { *triangle = malloc(sizeof(struct tr));}
        else
           *triangle = (struct tr*)realloc(*triangle, (*ntr + 1) * sizeof(struct tr));

        (*triangle + *ntr)->A = A; (*triangle + *ntr)->B = B; (*triangle + *ntr)->C = C;
        *ntr += 1; ++i1; ++i2;

        printf("Треугольник №%d добавлен успешно!\n", *ntr);
    }
}


void SPolyDiv(int *poly, int npoly, struct tr* triangle, int *ntr)
{
    int i0 = 0, i1 = 1, i2 = 2, A, B, C;

    while(i2 != npoly + 1) {
        A = poly[i0], B = poly[i1], C = poly[i2];

    //    RememberVertex(A, B, C, ntr);
        printf("%d\n", *ntr);
        if (*ntr == 0)
            { triangle = malloc(sizeof(struct tr));}
        else
            triangle = (struct tr*)realloc(triangle, (*ntr + 1) * sizeof(struct tr));

        (triangle + *ntr)->A = A; (triangle + *ntr)->B = B; (triangle + *ntr)->C = C;
        *ntr += 1; ++i1; ++i2;

        printf("Треугольник №%d добавлен успешно!\n", *ntr);
    }
}


void extra_init()
{

    int i;
    zbuffer = (int *)malloc((Xvp_range + 1)* (Yvp_range + 1)* sizeof(int));
    for (i = 0; i < (Xvp_range + 1) * (Yvp_range + 1); ++i)
            *(zbuffer + i) = 20000;
}

void ReadVertex(char *str, struct ver **fvertex, int *nvr)
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


    if (*nvr == 1) {
	*fvertex = (struct ver *)malloc(2 * sizeof(struct ver));
	(*fvertex)[0].x = 0; (*fvertex)[0].y = 0; (*fvertex)[0].z = 0;
    } else
        *fvertex = (struct ver *)realloc(*fvertex, (*nvr + 1) * sizeof(struct ver));

        (*fvertex + *nvr)->x = x; (*fvertex + *nvr)->y = y; (*fvertex + *nvr)->z = z;
        *nvr += 1;

}

/*void ReadNorm(char *str, int num)
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

}*/

void CalcVertex(struct ver* fvertex, struct ver** vertex, int nvr, double *Range)
{
    int i;
    double x, y, z, xe, ye, ze, X, Y;
    double xO = fvertex[0].x, yO = fvertex[0].y, zO = fvertex[0].z;
    for (i = 1; i <= nvr; ++i) {
        x = (fvertex + i)->x; y = (fvertex + i)->y; z = (fvertex + i)->z;

        //преобразуем координаты (через матрицу поворота)
        viewing(x,xO, y,yO, z,zO, &xe, &ye, &ze);
//        if (ze <= eps)
//            continue;

        //находим max и min для масштабирования в наше окно вывода
        X = xe / ze; Y = ye / ze;

        if (X < Range[0]) Range[0] = X; if (X > Range[1]) Range[1] = X;   //Range = {Xmin, Xmax, Ymin, Ymax}
        if (Y < Range[2]) Range[2] = Y; if (Y > Range[3]) Range[3] = Y;

        (*vertex + i)->x = xe; (*vertex + i)->y = ye; (*vertex + i)->z = ze;

    }
}

void extraCalcVertex(struct ver** vertex, int nvr, struct tr* triangle)
{
    int i, j, *ptr, Ntr;
    double a, b, c, r;
    for (i = 1; i <= nvr; ++i) {
	a = 0; b = 0; c = 0;
        ptr = (*vertex + i)->connect;

        Ntr = *ptr;

        for (j = 1; j <= Ntr; ++j) {
            a += (triangle + *(ptr + j) )->a;
            b += (triangle + *(ptr + j) )->b;
            c += (triangle + *(ptr + j) )->c;
        }

        r = sqrt(a*a + b*b + c*c);
        a /= r; b /= r; c /= r;

        (*vertex + i)->a = a; (*vertex + i)->b = b; (*vertex + i)->c = c;
    }

}

void CalcTr(struct ver* vertex, struct tr **triangle, int ntr)
{
    double xA, yA, zA, xB, yB, zB, xC, yC, zC, a, b, c, h, r;
    int i, A, B, C;
    for (i = 0; i < ntr; ++i) {
        A = (*triangle + i)->A; B = (*triangle + i)->B; C = (*triangle + i)->C;

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

        (*triangle + i)->a = a; (*triangle + i)->b = b; (*triangle + i)->c = c;
	(*triangle + i)->h = h;
    }
}

void scale_calc(double Xmin, double Xmax, double Ymin, double Ymax)
{
    double fx, fy, Xcentre, Ycentre, Xvp_centre, Yvp_centre;
    double Xrange, Yrange;

    Xrange = Xmax - Xmin; Yrange = Ymax - Ymin;
    Xvp_range = Xvp_max - Xvp_min; Yvp_range = Yvp_max - Yvp_min;
    fx = Xvp_range / Xrange; fy = Yvp_range / Yrange;
    d = (fx < fy) ? fx : fy;
    Xcentre = (Xmin + Xmax) / 2; Ycentre = (Ymax + Ymin) / 2;
    Xvp_centre = (Xvp_min + Xvp_max) / 2; Yvp_centre = (Yvp_min + Yvp_max) / 2;
    c1 = Xvp_centre - d * Xcentre; c2 = Yvp_centre - d * Ycentre;
}

void read_and_draw(struct ver** fvertex, struct ver** vertex, int nvr, struct tr** triangle, int ntr, double rho, double theta, double phi)
{

    //update();
    double Range[4] = {BIG, -BIG, BIG, -BIG};
    coeff(rho, theta, phi);

    CalcVertex(*fvertex, vertex, nvr, Range);

    if (onetime) {
        scale_calc(Range[0], Range[1], Range[2], Range[3]);
	extra_init();
    }
    onetime = 0;

    CalcTr(*vertex, triangle, ntr);

    /* Считаем нормали к вершинам*/
    extraCalcVertex(vertex, nvr, *triangle);

    DrawTr(*triangle, ntr, *vertex);
    DrawButtons();

}

void RememberVertex(struct tr* triangle, int ntr, struct ver** vertex)
{
    int *ptr, n, A, B, C, i;
    for (i = 0; i < ntr; ++i) {
        A = (triangle + i)->A; B = (triangle + i)->B; C = (triangle + i)->C;

        ptr = (*vertex + A)->connect; n = *ptr; ++n;
// ptr = (int *)realloc(ptr, (n + 1) * sizeof(int));
        *(ptr + n) = i; *ptr = n;

        ptr = (*vertex + B)->connect; n = *ptr; ++n;
// ptr = (int *)realloc(ptr, (n + 1) * sizeof(int));
        *(ptr + n) = i; *ptr = n;

        ptr = (*vertex + C)->connect; n = *ptr; ++n;
// ptr = (int *)realloc(ptr, (n + 1) * sizeof(int));
        *(ptr + n) = i; *ptr = n;
    }
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

void viewing(double x, double xO,double y, double yO, double z, double zO, double *pxe, double *pye, double *pze)
{
    *pxe = (v11 * (x - xO) + v21 * (y - yO)) + xO;
    *pye = (v12 * (x - xO) + v22 * (y - yO) + v32 * (z - zO)) + yO;
    *pze = (v13 * (x - xO)+ v23 * (y - yO) + v33 * (z - zO) + v43) + zO;

}

void FinalizeMe(struct ver** vertex, int nvr, struct ver** fvertex, struct tr** triangle)
{
    int i;

    free(zbuffer);
    free(*triangle);
    for (i = 1; i <= nvr; ++i)
        free((*vertex)[i].connect);
    free(*vertex);
    free(*fvertex);
}
