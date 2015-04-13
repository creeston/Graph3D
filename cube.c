#include<stdio.h>
#include<graphics.h>
#include<math.h>
#include<time.h>
//#include<X11/Xlib.h>
/*
* TODO:
* Матрицы в структуру и передавать в функцию
* посмотреть ошибку
* Удаление невидимых линий
*/
float v11, v12, v13,
      v21, v22, v23,
      v31, v32, v33, v43,
      screen_dist = 100, c1 = 220.0, c2 = 180.0, h = 250.0, rho = 1000, theta = 30, phi = 30;

int gd = DETECT, gm;

FILE *fp;

void coeff(float rho, float theta, float phi);
void mv(float x, float y, float z);
void dw(float x, float y, float z);
void perspective(float x, float y, float z, float *pX, float *pY);
void draw_from_file(char *filename);
int check(char ch);
void update(char *filename);

int main(int argc, char *argv[])
{
    printf("start? huh \n");
    char c = getchar();

    if (argc != 2) {
	printf("Invalid output! \n ./cube <filename> \n");
	return 0;
    }

    initgraph(&gd, &gm, NULL);

    coeff(rho, theta, phi);

    draw_from_file(argv[1]);

    while (1 == 1) {
        if ( check('w') )
	    for (; !kbhit(); phi -= 3.0)
                update(argv[1]);
        if (check('s') )
	    for (; !kbhit(); phi += 3.0)
	        update(argv[1]);
        if (check('a') )
            for (; !kbhit(); theta += 3.0)
                update(argv[1]);
        if (check('d') )
            for (; !kbhit(); theta -= 3.0)
                update(argv[1]);
        if (check('q') )
            for (; !kbhit(); screen_dist -= 4.0)
                update(argv[1]);
        if (check('e') )
            for (; !kbhit(); screen_dist += 4.0)
                update(argv[1]);
        if (check('x') )
            closegraph();
    }
}


void coeff(float rho, float theta, float phi)
{
    float th, ph, costh, sinth, cosph, sinph, factor;

    factor = atan(1.0) / 45.0;
    th = theta * factor; ph = phi * factor;
    costh = cos(th); sinth = sin(th);
    cosph = cos(ph); sinph = sin(ph);

    v11 = -sinth; v12 = -cosph * costh; v13 = -sinph * costh;
    v21 = costh;  v22 = -cosph * sinth; v23 = -sinph * sinth;
    v31 = 0;      v32 = sinph;          v33 = -cosph;        v43 = rho;
}


void mv(float x, float y, float z)
{
    float X, Y;
    perspective(x,y,z, &X, &Y);
    moveto((int)X,(int)Y);
}

void dw(float x, float y, float z)
{
    float X, Y;
    perspective(x,y,z, &X,  &Y);
    lineto((int)X, (int)Y);
}

void perspective(float x, float y, float z, float *pX, float *pY)
{
    float xe, ye, ze;
    xe = v11 * x + v21 * y;
    ye = v12 * x + v22 * y + v32 * z;
    ze = v13 * x + v23 * y + v33 * z + v43;

    *pX = screen_dist * xe / ze + c1;
    *pY = screen_dist * ye / ze + c2;
}


int check(char ch)
{
    char c;
    if (kbhit()) c = getchar();
    if (c == ch) return 1;
    else return 0;
}

void update(char *filename)
{
    cleardevice();
    coeff(rho, theta, phi);
    draw_from_file(filename);
    delay(1);
}

void draw_from_file(char *filename)
{
    fp = fopen(filename, "r");
    if (fp == NULL && filename != NULL) {
	printf("This file doesn't exist");
	exit(0);
    }

    float x, y, z;
    int i;

    while ( fscanf(fp, "%f %f %f %d", &x,&y,&z,&i) > 0) {
	if (i)
            dw(x,y,z);
        else mv(x,y,z);
    }
    fclose(fp);
}

