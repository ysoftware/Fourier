raylib_i = $(shell pkg-config --cflags raylib)
raylib_l = $(shell pkg-config --libs raylib)

all:
	cc -o main $(raylib_i) main.c -lm $(raylib_l)
	./main
