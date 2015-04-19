#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NVERTEX 50
#define NTRIANGLE 50
/*
*
*   TODO
* 1. make it works
* 2. link with main.c
* 3. make it rational
*/
int ntr = 0;
double v11, v12, v13, v21, v22, v23, v32, v33, v43,
	eps = 1e-5, meps = -1e-5, oneminus = 1 - 1.e-5, oneplus = 1 + 1.e-5;

FILE *fpout;
struct {
    float X,Y;
    int code;
} s;

struct {
    double x,y,z;
} VERTEX[NVERTEX];

struct {
    int A,B,C;
    double a,b,c,h;
} TRIANGLE[NTRIANGLE];

int error(char *str)
{
    printf("%s \n", str);
    exit(1);
}

void check()
{
    static int i = 1;
    printf("input to continue");
    char c = getchar();
    printf("Checkpoint # %d", i);
    i++;
}
int main(int argc, char *argv[])
{
    int A,B,C,P,Q;
    double xO, yO, zO, rho, theta, phi, x, y, z, xA, yA, zA, xB, yB, zB, xC, yC, zC, a, b, c, h, r;
    FILE *fpin;
    int ch, i;

    if (argc != 2 || (fpin = fopen(argv[1], "r")) == NULL) {
        error("Incorrect input file \n");
    }

    fscanf(fpin, "%lf %lf %lf", &xO,&yO,&zO); skipf(fpin);
    fscanf(fpin, "%lf %lf %lf", &rho,&theta,&phi);
    coeff(rho, theta, phi);

    check(); //1

    while(ch = getc(fpin), ch != '#' && ch != EOF) {
        ungetc(ch, fpin);
        fscanf(fpin, "%d %lf %lf %lf", &i, &x, &y, &z);
        if (i < 0 || i > NVERTEX) {
            error("Invalid vertex number \n");
        }
        viewing(x - xO, y - yO, z - zO, &VERTEX[i].x, &VERTEX[i].y, &VERTEX[i].z);
        if (VERTEX[i].z < eps)
            error("deifferent sides of viepoint \n");

    }
    check(); //2

    while(ch = getc(fpin), ch != '#') {
        ungetc(ch, fpin);
        fscanf(fpin, "%d %d %d", &A,&B,&C);
        xA = VERTEX[A].x; yA = VERTEX[A].y; zA = VERTEX[A].z;
        xB = VERTEX[B].x; yB = VERTEX[B].y; zB = VERTEX[B].z;
        xC = VERTEX[C].x; yC = VERTEX[C].y; zC = VERTEX[C].z;

        a = yA * (zB - zC) - yB * (zA - zC) + yC * (zA - zB);
        b = -(xA * (zB - zC) - xB * (zA - zC) + xC * (zA - zB));
        c =  xA * (yB - yC) - xB * (yA - yC) + xC * (yA - yB);
        h = xA * (yB * zC - yC * zB) -
            xB * (yA * zC - yC * zA) +
            xC * (yA * zB - yB * zA);
        if (h > 0) {
            if (ntr == NTRIANGLE) {
		printf("%d \n", ntr);
                error("Too much triangles");
	    }

            r = sqrt(a * a + b * b + c * c);
            a /= r; b /= r; c /= r;
            TRIANGLE[ntr].A = A;
            TRIANGLE[ntr].B = B;
            TRIANGLE[ntr].C = C;
            TRIANGLE[ntr].a = a;
            TRIANGLE[ntr].b = b;
            TRIANGLE[ntr].c = c;
            TRIANGLE[ntr++].h = h;
        }
    }
    check(); //3

    fpout = fopen("a.scratch", "w");
    if (fpout == NULL)
	error("File cannot be opened");
    while(fscanf(fpin, "%d %d", &P, &Q) > 0) {
        linesegment(VERTEX[P].x, VERTEX[P].y, VERTEX[P].z, VERTEX[Q].x, VERTEX[Q].y, VERTEX[Q].z, 0);
    }
    check(); //4
    fclose(fpout);
}

skipf(FILE *fpin)
{
    int ch;
    while (ch = getc(fpin), ch != '\n' && ch !=  EOF);
}

coeff(double rho, double theta, double phi)
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

viewing(double x, double y, double z, double *pxe, double *pye, double *pze)
{
    *pxe = v11 * x + v21 * y;
    *pye = v12 * x + v22 * y + v32 * z;
    *pze = v13 * x + v23 * y + v33 * z + v43;
}

linesegment(double xP, double yP, double zP, double xQ, double yQ, double zQ, int j0)
{
    int j = j0, worktodo = 1, A, B, C, i, Pbeyond, Qbeyond, outside, Poutside, Qoutside, eA, eB, eC, sum;
    double xA, yA, zA, xB, yB, zB, xC, yC, zC,
	   dA, dB, dC, MIN, MAX, lab, mu,
	   xmin, ymin, zmin, xmax, ymax, zmax,
	   C1, C2, C3, K1, K2, K3, denom1, denom2,
	   Cpos, Ppos, Qpos, aux, eps1, a, b, c, h, hP, hQ, r1, r2, r3;

    while(j < ntr) {
        a = TRIANGLE[j].a; b = TRIANGLE[j].b;  c = TRIANGLE[j].c;
        h = TRIANGLE[j].h; eps1 = eps + eps * h;


       /* TEST #1*/

        hP = a * xP + b * yP + c * zP;
        hQ = a * xQ + b * yQ + c * zQ;
        if (hP < h + eps1 && hQ < h + eps1) {
            ++j;          /* PQ is not behind ABC*/
            continue;
        }

        /* TEST #2*/

        K1 = yP * zQ - yQ * zP; K2 = zP * xQ - zQ * xP; K3 = xP * yQ - xQ * yP;
        A = TRIANGLE[j].A;  B = TRIANGLE[j].B; C = TRIANGLE[j].C;
        xA = VERTEX[A].x; yA = VERTEX[A].y; zA = VERTEX[A].z;
        xB = VERTEX[B].x; yB = VERTEX[B].y; zB = VERTEX[B].z;
        xC = VERTEX[C].x; yC = VERTEX[C].y; zC = VERTEX[C].z;

        dA = K1 * xA + K2 * yA + K3 * zA;
        dB = K1 * xB + K2 * yB + K3 * zB;
        dC = K1 * xC + K2 * yC + K3 * zC;

        eA = dA > eps ? 1 : dA < meps ? -1 : 0;
        eA = dB > eps ? 1 : dB < meps ? -1 : 0;
        eA = dC > eps ? 1 : dC < meps ? -1 : 0;

        sum = eA + eB + eC;
        if (abs(sum) >= 2) {
            ++j;
            continue;
        }

        /* TEST #3*/

        Poutside = Qoutside = 0; MIN = 1; MAX = 0;
        for(i = 0; i < 3; ++i) {
            C1 = yA * zB - yB * zA; C2 = zA * zB - zB * xA; C3 = xA * yB - xB * yA;
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
            lab = fabs(denom1) <= eps ? 1.e7 : -Ppos / denom1;

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
            j++;
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
            ++j;
            continue;
        }

        xmax = xP + MAX * r2; ymax = yP + MAX * r2; zmax = zP + MAX * r3;
        if (a * xmax + b * ymax + c * zmax < h - eps1) {
            ++j;
            continue;
        }

        /* TEST #6*/

        if (Poutside)
            linesegment(xP, yP, zP, xmin, ymin, zmin, j + 1);
        if (Qoutside)
            linesegment(xQ, yQ, zQ, xmax, ymax, zmax, j + 1);
        worktodo = 0;
        break;
    }

    /*if (worktodo) {
        s.X = xP / zP; s.Y = yP / zP; s.code = 0;
        fwrite(&s, sizeof s, 1, fpout);
        s.X = xQ / zQ; s.Y = yQ / zQ; s.code = 1;
        fwrite(&s, sizeof s, 1, fpout);
    }*/

    if (worktodo) {
        fprintf(fpout, "%f %f %d\n", xP / zP * 400, yP / zP * 400, 0);
	fprintf(fpout, "%f %f %d\n", xQ / zQ * 400, yQ / zQ * 400, 1);
    }
}


