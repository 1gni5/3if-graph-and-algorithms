# simulated annealing applied to n-queens

# CFLAGS = -O3  (for best performance)
# CFLAGS = -g   (for debugging)
CFLAGS = -O3 -Wall

LIBS = -lm

tsp: tsp.o sa.o
	cc $(CFLAGS) -o tsp tsp.o sa.o $(LIBS)

tsp.o: tsp.c sa.h
	cc $(CFLAGS) -c tsp.c

nq: nq.o sa.o
	cc $(CFLAGS) -o nq nq.o sa.o $(LIBS)

nq.o: nq.c sa.h
	cc $(CFLAGS) -c nq.c

sa.o: sa.c sa.h
	cc $(CFLAGS) -c sa.c

clean:
	rm -f nq.o sa.o nq
