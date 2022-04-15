#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* https://stackoverflow.com/questions/1644868/ */
#define debug(fmt, ...) \
          do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)

#define DEBUG 1

void scanParam();
int scanConfigTitle();
void allocTempsAndCosts();
void sa();
void heat();
void initPlot();
void plot();

int cost();
void rndMove();
void undoMove();
void saveNewBest();
