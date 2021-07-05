CC      = gcc-9
CFLAGS	= -Ofast -march=native -ffast-math -funroll-loops -floop-block -ftree-vectorize #-fopt-info-vec-missed -fopt-info-vec
CPPFLAGS=
WFLAGS	= -std=gnu11 -Wall -Wextra -g
LDFLAGS	= -lm

ISPC = /opt/ispc/1.15.0/bin/ispc
ISPCFLAGS = -O3 --math-lib=fast --target avx2-i64x4 #fast-masked-vload fast-math


# TARGETS	= tiny_md viz
TARGETS	= tiny_md
SOURCES	= $(shell echo *.c)
OBJECTS = vec.o core.o

all: $(TARGETS)

# viz: viz.o $(OBJECTS)
#	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lGL -lGLU -lglut


tiny_md: tiny_md.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

vec.o: vec.ispc
	$(ISPC) $(ISPCFLAGS) $(CPPFLAGS) -h vec.h -o $@ $^

%.o: %.c 
	$(CC) $(WFLAGS) $(CPPFLAGS) $(CFLAGS) -c $<


clean:
	rm -f $(TARGETS) *.o *.xyz *.log .depend *.temp 

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend

.PHONY: clean all

