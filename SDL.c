#include <stdio.h>
#include <SDL/SDL.h>
#include <stdlib.h>
#include <unistd.h>

#include "test_header.h"

SDL_Surface *screen;

void clearscreen()
{
    int i,j;
    struct color Cl = {255,255,255};
    SDL_UpdateRect(screen,0,0,0,0);
    SDL_FillRect(screen, NULL, 0x000000);
}

void init(int width, int height, int depth)
{

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("Couldn't initialize SDL: %s \n", SDL_GetError());
        exit(1);
    }

    atexit(SDL_Quit);

    //SDL_DisplayFormat(screen);
    screen = SDL_SetVideoMode(width, height, depth, SDL_SWSURFACE | SDL_ANYFORMAT);
    if (screen == NULL)
	printf("error\n");

}

void drawpix(int x, int y, struct color C)
{
    Uint32 dye;
    dye = SDL_MapRGB(screen->format, C.r, C.g, C.b);

/*    if (SDL_MUSTLOCK(screen)) {
        if (SDL_LockSurface(screen) < 0) {
            fprintf(stderr, "Cant lock screen: %s \n", SDL_GetError());
            return;
        }
    }*/

    putpix(x, y, dye);

  /*  if (SDL_MUSTLOCK(screen)) {
        SDL_UnlockSurface(screen);
    }*/

    //SDL_UpdateRect(screen, x, y, 1, 1);

}

/* int main(int argc, char *argv[])
{
    int width, height, x;

    scanf("%d %d", &width, &height);
    init(width, height, 8);
    struct color Cl = {217,255,147};
    for (x = 1; x <= width; ++x) {
        drawpix(x, height / 2, Cl);
	sleep(1);
    }
    return 0;
}*/

void putpix(int x, int y, Uint32 pixel)
{
    int bpp = screen->format->BytesPerPixel;
    Uint8 *p = (Uint8 *)screen->pixels + y * screen->pitch + x * bpp;
    //p = pixel;
    switch (bpp) {
        case 1:
            *p = pixel;
            break;
        case 2:
            *(Uint16 *)p = pixel;
            break;
        case 3:
            if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
                p[0] = (pixel >> 16) & 0xff;
                p[1] = (pixel >> 8) & 0xff;
                p[2] = pixel & 0xff;
            } else {
                p[0] = pixel & 0xff;
                p[1] = (pixel >> 8) & 0xff;
                p[2] = (pixel >> 16) & 0xff;
            }
            break;
        case 4:
            *(Uint32 *)p = pixel;
            break;
    }
    
}

