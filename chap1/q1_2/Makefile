.PHONY: clean, all, show

CC      := gcc
CFLAGS  := -O2 -std=c11 -Wall
LDFLAGS := -L ../../lib/ -Wl,-rpath=../../lib/ -lm -lmt

ifdef NAME
all: $(NAME).png
else
all: a.png
endif

%.out: %.c Makefile
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

%.dat: %.out Makefile
	./$< > $(*F).dat

%.png: %.plt %.dat Makefile
	gnuplot $<

ifdef NAME
show: $(NAME).png Makefile
else
show: a.png Makefile
endif
	gnome-open $<

clean:
	-rm -f *.o *.out *.dat *.png
