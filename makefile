all:
	gcc hidlinpx.c linesegment.c -g -lm -lgraph -lX11 -lSDL -o PROGRAMM
