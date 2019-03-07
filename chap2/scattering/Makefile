.PHONY: show

CPPFLAGS := -O2 -std=gnu++17 -Wall -fopenmp
LDFLAGS  := -L ../../lib/ -Wl,-rpath=../../lib/ -lode -lgsl -lblas

cross_section.out: cross_section.cpp Makefile
	g++ $(CPPFLAGS) $< -o $@ $(LDFLAGS)

cross_section.dat: cross_section.out Makefile
	./$< > $@

cross_section.png: cross_section.plt cross_section.dat Makefile
	gnuplot $<

show: cross_section.png Makefile
	gnome-open $<