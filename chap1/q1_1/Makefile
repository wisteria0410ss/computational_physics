.PHONY: clean, all, show

CC      := gcc
CFLAGS  := -O2 -std=c11 -Wall
LDFLAGS := -lm
OBJS    := rk4.o duffing.o

ifdef NAME
all: $(NAME).png
else
all: a.png
endif

%.o: %.c duffing.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

%.out: %.c $(OBJS) duffing.h Makefile
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

%.dat: %.out Makefile
	./$< 0.0 1.0 > $(*F).dat
	cat $(*F).dat

%.png: %.plt %.dat Makefile
	gnuplot $<

ifdef NAME
show: $(NAME).png Makefile
else
show: a.png Makefile
endif
	gnome-open $<

clean:
	-rm *.o *.out *.dat *.png