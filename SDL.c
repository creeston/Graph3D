/*
*  TODO
*  SDL -> SDL2
*     Проблемы, из которых я ещё этого не сделал: 
*	1. Не воспринимает SDL_Eventы в main функции
*
*/
#include <stdio.h>
#include <SDL/SDL.h>
#include <stdlib.h>
#include <unistd.h>

#include "test_header.h"

//SDL_Window *window;
//SDL_Renderer *render;
SDL_Surface *screen;
SDL_Surface *pic;
void clearscreen()
{
    //SDL_SetRenderDrawColor(render, 0,0,0,255);
    //SDL_RenderClear(render);
    //SDL_RenderPresent(render);
    SDL_UpdateRect(screen,0,0,0,0);
    SDL_FillRect(screen, NULL, 0x000000);
}

void update_render()
{
//    SDL_RenderPresent(render);
}

void init(int width, int height, int depth)
{
//    SDL_Init(SDL_INIT_EVERYTHING);
//    window = SDL_CreateWindow("3d", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 640, 480, SDL_WINDOW_SHOWN);
//    render = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
//    SDL_SetRenderDrawColor(render, 0,0,0,255);
//    SDL_RenderClear(render);
//    SDL_RenderPresent(render);

    screen = SDL_SetVideoMode(width, height, depth, SDL_SWSURFACE | SDL_ANYFORMAT);
    pic = load_image("grass2.bmp");
}

void drawpix(int x, int y, struct color C)
{
    Uint32 dye;
    dye = SDL_MapRGB(screen->format, C.r, C.g, C.b);
    //put_pixel32(x, y, dye);
    putpix(x, y, dye);
}

SDL_Event GetEvent()
{
    SDL_Event event;
    SDL_PollEvent(&event);
    return event;
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
void LoadTexture(SDL_Surface *image, int x, int y)
{
    //SDL_Surface *image = load_image(filename);
    Uint32 pixel = get_pixel32(pic, x ,y);
    put_pixel32(x, y, pixel);
    //SDL_FreeSurface(image);
}

SDL_Surface *load_image(char *str)
{
    SDL_Surface *loadedImage = NULL;
    SDL_Surface *optimizedImage = NULL;

    loadedImage = SDL_LoadBMP(str);
    if (loadedImage != NULL) {
        optimizedImage = SDL_DisplayFormat(loadedImage);
        SDL_FreeSurface(loadedImage);
    }

    return optimizedImage;
}

void applySurface(int x, int y, int w, int h, SDL_Surface *source)
{
    SDL_Rect offset;
    offset.x = x;
    offset.y = y;
    offset.h = h;
    offset.w = w;

    SDL_BlitSurface(source, NULL, screen, &offset);
}


Uint32 get_pixel32(SDL_Surface *surface, int x, int y)
{
    Uint32 *pixels = (Uint32 *)surface->pixels;
    return pixels[y * surface->w + x];
}

void put_pixel32(int x, int y, Uint32 pixel)
{
    Uint32 *pixels = (Uint32 *)screen->pixels;
    pixels[y * screen->w + x] = pixel;
}

void line(int x0, int y0, int x1, int y1, struct color Cl)
{
    int steep = 0;
    if (abs(x0 - x1) < abs(y0 - y1)) {
        swapI(&x0, &y0);
        swapI(&x1, &y1);
        steep = 1;
    }

    if (x0 > x1) {
        swapI(&x0, &x1);
        swapI(&y0, &y1);
    }

    for (int x = x0; x <= x1; ++x) {
        float t = (x - x0) / (float)(x1 - x0);
        int y = y0 * (1. - t) + y1 * t;
        if (steep)
            drawpix(y, x, Cl);
        else
            drawpix(x, y, Cl);
    }
}

void DrawRect(int x1, int y1, int x2, int y2, struct color C, int width)
{
    for (int i = 0; i < width; ++i) {
        line(x1 - i, y1 - i, x2 + i, y1 - i, C);
        line(x2 + i, y1 - i, x2 + i, y2 + i, C);
        line(x2 + i, y2 + i, x1 - i, y2 + i, C);
        line(x1 - i, y2 + i, x1 - i, y1 + i, C);
    }
}
