#include<stdio.h>
#include<graphics.h>
#include<math.h>
#include<time.h>
//#include<X11/Xlib.h>
/*
* TODO:
* посмотреть ошибку
* Удаление невидимых линий
*/
float screen_dist = 100, c1 = 220.0, c2 = 180.0, h = 250.0, rho = 1000, theta = 30, phi = 30;

float **matrix;

FILE *fp;
struct vertex {
    float x;
    float y;
    float z;
    int i;
};

struct vertex **obj;

void opengraph();
void coeff(float rho, float theta, float phi);
void mv(float x, float y, float z);
void dw(float x, float y, float z);
void perspective(float *cords, float *pX, float *pY);
void draw_from_file(char *filename);
int check(char ch);
void update(int n);
void mx_mult(float *mx1, float **mx2, float *result, int a, int b);
void checkp();
void copycord(char *filename, int *m);
void draw_from_object(int n);
void merge(char *filename, int *n);
void checkp()
{
    static int count = 1;
    printf("Checkpoint # %d", count);
    getchar();
    ++count;
}
void mx_mult(float *mx1, float **mx2, float *result, int a, int b)
{
    int i, j;
    float sum;;
    for(i = 0; i <= b - 1; ++i) {
        for (j = 0, sum = 0; j <= a - 1; ++j) {
            sum += mx1[j] * mx2[j][i];
        }
        result[i] = sum;
    }
}

int main(int argc, char *argv[])
{
    int i, num;

    matrix = malloc(4 * sizeof(float *));
    for (i = 0; i < 4; ++i) {
	matrix[i] = malloc(4 * sizeof(float));
    }

    printf("start? huh \n");
    if (argc < 2) {
	printf("Invalid output! \n ./cube <filename> \n");
	return 0;
    }

    coeff(rho, theta, phi);
    copycord(argv[1], &num);

    opengraph();
    draw_from_object(num);

    while (1 == 1) {
        if ( check('w') )
	    for (; !kbhit(); phi -= 3.0)
                update(num);
        if (check('s') )
	    for (; !kbhit(); phi += 3.0)
	        update(num);
        if (check('a') )
            for (; !kbhit(); theta += 3.0)
                update(num);
        if (check('d') )
            for (; !kbhit(); theta -= 3.0)
                update(num);
        if (check('q') )
            for (; !kbhit(); screen_dist -= 4.0)
                update(num);
        if (check('e') )
            for (; !kbhit(); screen_dist += 4.0)
                update(num);
	if (check('m') )
	    merge(argv[2], &num);
        if (check('x') )
            goto exit;
    }
    exit:
    closegraph();
    /*Free memory */
    for (i = 0; i < num; ++i) {
	free(obj[i]);
    }
    free(obj);
    for (i = 0; i < 4; ++i) {
	free(matrix[i]);
    }
    free(matrix);
}

void opengraph()
{
    int gd = DETECT, gm;
    initgraph(&gd, &gm, NULL);
}

void coeff(float rho, float theta, float phi)
{
    float th, ph, costh, sinth, cosph, sinph, factor;

    factor = atan(1.0) / 45.0;
    th = theta * factor; ph = phi * factor;
    costh = cos(th); sinth = sin(th);
    cosph = cos(ph); sinph = sin(ph);

    matrix[0][0] = -sinth; matrix[0][1] = -cosph * costh; matrix[0][2] = -sinph * costh; matrix[0][3] = 0.0;
    matrix[1][0] = costh;  matrix[1][1] = -cosph * sinth; matrix[1][2] = -sinph * sinth; matrix[1][3] = 0.0;
    matrix[2][0] = 0;      matrix[2][1] = sinph;          matrix[2][2] = -cosph;         matrix[2][3] = 0.0;
    matrix[3][0] = 0;      matrix[3][1] = 0;              matrix[3][2] = rho;            matrix[3][3] = 1.0;
 }


void mv(float x, float y, float z)
{
    float X, Y;
    float cords[] = {x, y, z, 1.0};
    perspective(cords, &X, &Y);
    moveto((int)X,(int)Y);
}

void dw(float x, float y, float z)
{
    float X, Y;
    float cords[] = {x, y, z, 1.0};
    perspective(cords, &X,  &Y);
    lineto((int)X, (int)Y);
}

void perspective(float *cords, float *pX, float *pY)
{
    float new_cords[4];
    mx_mult(cords, matrix, new_cords, 4, 4);
    *pX = screen_dist * new_cords[0] / new_cords[2] + c1;
    *pY = screen_dist * new_cords[1] / new_cords[2] + c2;
}


int check(char ch)
{
    char c;
    if (kbhit()) c = getchar();
    if (c == ch) return 1;
    else return 0;
}

void update(int n)
{
    cleardevice();
    coeff(rho, theta, phi);
    draw_from_object(n);
    delay(1);
}

void copycord(char *filename, int *m)
{
    fp = fopen(filename, "r");
    if (fp == NULL && filename != NULL) {
        printf("This file doesn't exist");
        exit(0);
    }
    float x, y, z;
    int i,j = 0;
    (*m) = 0;

    while (fscanf(fp, "%f %f %f %d", &x,&y,&z,&i) > 0)
	(*m)++;
    fclose(fp);
    free(obj); //if memory was already allocated
    obj = malloc((*m) * sizeof(struct vertex *));
    for (j = 0; j < (*m); ++j) {
	obj[j] = malloc(sizeof(struct vertex));
    }

    fp = fopen(filename, "r");

    for (j = 0; fscanf(fp, "%f %f %f %d", &x,&y,&z,&i) > 0; ++j) {
        obj[j]->x = x;
	obj[j]->y = y;
	obj[j]->z = z;
	obj[j]->i = i;
    }
    fclose(fp);
}

void draw_from_object(int n)
{
    int j;
    for (j = 0; j < n; ++j) {
        if (obj[j]->i == 1)
            dw(obj[j]->x, obj[j]->y, obj[j]->z);
        else mv(obj[j]->x, obj[j]->y, obj[j]->z);
    }

}

void merge(char *filename, int *n)
{
    int i, j, n1;
    struct vertex **obj1;

    obj1 = malloc(*n * sizeof(struct vertex *));
    for (i = 0; i < *n; ++i)
	obj1[i] = malloc(sizeof(struct vertex));

    for (i = 0; i < *n; ++i)
	*obj1[i] = *obj[i];


    copycord(filename, &n1);


    obj1 = realloc(obj1, n1 * sizeof(struct vertex *));

    for (i = *n; i < *n + n1; ++i) {
        obj1[i] = malloc(sizeof(struct vertex));
    }

    for (i = *n, j = 0; i < *n + n1; ++i, ++j)
	*obj1[i] = *obj[j];

    free(obj);
    obj = obj1;
    *n += n1;
}
