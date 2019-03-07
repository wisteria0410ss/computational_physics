.PHONY: clean, all, show

CC      := gcc
CFLAGS  := -O2 -std=c11 -Wall
LDFLAGS := -L ../../lib/ -Wl,-rpath=../../lib/ -lm -lmt

all: a.png

a.o: a.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

a.out: a.c Makefile
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

a.dat: a.out Makefile
	./$< > a.dat

a.png: a.plt a.dat Makefile
	gnuplot $<

show: a.png Makefile
	gnome-open $<

clean:
	-rm -f *.o *.out *.dat *.png
