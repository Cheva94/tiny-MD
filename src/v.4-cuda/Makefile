CU       = nvcc
# CUFLAGS  = -O3 -Xcompiler=-Wall -Xcompiler=-Wextra -Xcompiler=-Ofast -Xcompiler=-march=native -Xcompiler=-ffast-math -Xcompiler=-funroll-loops -Xcompiler=-floop-block -Xcompiler=-ftree-vectorize -arch=sm_75 #-Xcompiler=-g
CUFLAGS  = -O3 -Xcompiler=-Wall -Xcompiler=-Wextra -arch=sm_75 #-g -G #-Xcompiler=-rdynamic
LDFLAGS	 = -lm
CPPFLAGS =
# El device query dice que en jupiterace es 7.5 y zx81 es 8.6

PROG = Qdah
SRC	 = $(shell echo *.cu)
OBJ  = kernel.o core.o tiny_md.o

all: $(PROG)

Qdah: $(OBJ)
	$(CU) $(CUFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cu
	$(CU) $(CUFLAGS) $(CPPFLAGS) -o $@ -c $<

clean:
	rm -f $(PROG) *.o *.xyz *.log .depend *.temp

.depend: $(SRC)
	$(CU) -M $^ > $@

-include .depend

.PHONY: clean all
