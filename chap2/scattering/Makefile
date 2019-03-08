.PHONY: clean, show
.PRECIOUS: %.out %.dat

CPPFLAGS := -Ofast -std=gnu++17 -Wall -fopenmp
LDFLAGS  := -L ../../lib/ -Wl,-rpath=../../lib/ -lode -lgsl -lblas

%.out: %.cpp Makefile
	g++ $(CPPFLAGS) $< -o $@ $(LDFLAGS)

%.dat: %.out Makefile
	./$< > $@

%.png: %.plt %.dat Makefile
	gnuplot $<

show: cross_section.png Makefile
	gnome-open $<
	
show_d: differential_cross_section.png Makefile
	gnome-open $<

clean:
	-rm -f *.out *.dat *.png