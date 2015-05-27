all:
	gcc hidlinpx.c SDL.c -g -std=c99 -lm -lgraph -lX11 -lSDL -o PROGRAMM
