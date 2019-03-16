.PHONY: all

CFLAGS  := -O2 -Wall -Wextra
LDFLAGS := -L ../lib/ -Wl,-rpath=../lib/ -lm -lcmatrix -lblas -llapack

all: $(patsubst %.c,%.out,$(wildcard *.c))

%.out: %.c Makefile
	gcc $(CFLAGS) $< -o $@ $(LDFLAGS)
