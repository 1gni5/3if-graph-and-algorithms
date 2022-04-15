#include "sa.h"

static float T0;         /* initial temperature                             */
static float ALPHA;      /* T[n+1] = ALPHA * T[n] ; ALPHA < 1               */
static int MAXTIME;      /* annealing stops after at most MAXTIME           */
static float BETA;       /* M[n]= time spent at temperature T[n]            */
                         /* M[n+1]= BETA * M[n] ; BETA > 1                  */
static float BETA0;      /* M[0] = BETA0 * MAXTIME                          */
static int MAXCOST;      /* max of min of costs at successive temperatures  */
static int MINCOST;      /* min of all costs encountered                    */
static char TITLE[100];  /* title of the plot                               */
static float *TEMPS;     /* succesive temperature values to plot            */
static int *COSTS;       /* succesive best costs at each temperature        */
static int NBITER;       /* nb of iterations through the outer loop         */
static int REJECTED=0;   /* was the last saStep rejected?                   */

void /* the first lines of stdin must contain parameters */
scanParam() {
  if (scanf("ALPHA = %f\n", &ALPHA) != 1) {
     fprintf(stderr, "error scan ALPHA\n");
     exit(1);
  }
  if (scanf("BETA = %f\n", &BETA) != 1) {
    fprintf(stderr, "error scan BETA\n");
    exit(1);
  }
  if (scanf("BETA0 = %f\n", &BETA0) != 1) {
    fprintf(stderr, "error scan BETA0\n");
    exit(1);
  }
  if (scanf("MAXTIME = %i\n", &MAXTIME) != 1) {
    fprintf(stderr, "error scan MAXTIME\n");
    exit(1);
  }
}

int
scanConfigTitle() {
  if (fgets(TITLE, 100, stdin) == NULL) {
    fprintf(stderr, "error reading title configuration\n"); 
    exit(1);
  }
  strtok(TITLE, "\n");
  if (strcmp(TITLE, "END") == 0) return 0;
  return 1;
}

void
allocTempsAndCosts() {
  free(TEMPS); TEMPS = (float *)malloc(MAXTIME*sizeof(float));
  free(COSTS); COSTS = (int *)malloc(MAXTIME*sizeof(int));
}

double /* random number between 0 and 1 */
rand01() { return rand()/(double)RAND_MAX; }

int /* Simulated annealing. Evaluation of a random neighbor. */
saStep(int curCost, float temp) {
  int newCost, deltaCost;
  rndMove(); newCost = cost(); deltaCost = newCost - curCost; 
  if (deltaCost < 0) {
    curCost = newCost; REJECTED = 0;
    if (curCost < MINCOST) { MINCOST = curCost; saveNewBest(); }
  } else {
    if (rand01() < exp(- deltaCost / temp)) {
      curCost = newCost; REJECTED = 0;
    } else {
      undoMove(); REJECTED = 1;
    }
  }
  return curCost;
}

void /* simulated annealing */
sa() {
  int solutionFound=0;
  int curCost = cost();
  int bestCostAtTemp = MINCOST = curCost; saveNewBest();
  float temp = T0; 
  int elapsed = 0;
  int m = (int) floor(BETA0*MAXTIME); /* aMount of time spent at current temp */
  int timer = m;
  int nbReject;
  MAXCOST = 0; NBITER = 0;
  while (elapsed<MAXTIME && !solutionFound && temp>=0.1) {
    bestCostAtTemp = curCost;
    nbReject = 0;
    while (timer!=0) {
      curCost = saStep(curCost, temp); nbReject += REJECTED;
      if(curCost == 0) {bestCostAtTemp = curCost; solutionFound=1; break;}
      else if(curCost < bestCostAtTemp) bestCostAtTemp = curCost;
      timer -= 1;
    }
    elapsed += m;
    debug("At T=%f, elapsed: %i / %i, best: %i, rejections: %i\n", 
          temp, elapsed, MAXTIME, bestCostAtTemp, nbReject);
    TEMPS[NBITER] = temp; COSTS[NBITER] = bestCostAtTemp; NBITER += 1;
    m = (int) floor(m*BETA); /* More time spent at lower temp */
    timer = m; temp *= ALPHA;
    if(bestCostAtTemp > MAXCOST) MAXCOST = bestCostAtTemp;
  }
}

void      /* set T0 by starting at a low temperature and gradually heating */
heat() {  /* until the % of rejections is less than rejectionThreshold     */
  int curCost = cost();
  float temp = 1.0;
  float prctReject = 1.0;
  float rejectionThreshold = 0.05;
  int timer = 100;
  int nbReject;
  int i,k;
  int maxIter = 1000;

  for (k=0;k<maxIter && prctReject>rejectionThreshold;k++) {
    nbReject = 0;
    for (i=0;i<timer;i++) {
      curCost = saStep(curCost, temp);
      nbReject += REJECTED;
    }
    prctReject = (float)nbReject / (float)timer;
    debug("prct rejected at %f = %f\n", temp, prctReject);
    temp *= 1.1;
  }
  if (k>=maxIter) debug("heat() reaches maxIter at temp %f\n", temp);
  T0 = temp;
}

void /* copy the content of file s in file d */
fcopy(FILE *s, FILE *d) {int c; while ((c=getc(s))!=EOF) putc(c,d);}

void /* load a prelude of postscript functions */
initPlot() { 
  FILE *fpPlot = fopen("plot.ps", "w");
  if (fpPlot == NULL) {
    fprintf(stderr, "error fopen plot.ps in initPlot()\n");
    exit(1);
  }
  FILE *fpBase = fopen("base.ps", "r"); 
  if (fpBase == NULL) {
    fprintf(stderr, "error fopen base.ps in initPlot()\n");
    exit(1);
  }
  fcopy(fpBase, fpPlot); 
  fclose(fpBase); fclose(fpPlot);
}

void /* plot cost against temperature */
plot() {
  FILE *fpPlot = fopen("plot.ps", "a");
  if (fpPlot == NULL) {
    fprintf(stderr, "error fopen plot.ps in plot()\n");
    exit(1);
  }
  int i;
  int xmin = (int) floor(TEMPS[NBITER-1]), ymin = MINCOST;
  int xmax = (int) ceil(T0),               ymax = MAXCOST;
  fprintf(fpPlot, "newpath\n");
  fprintf(fpPlot, "/xmin %i def\n", xmin);
  fprintf(fpPlot, "/ymin %i def\n", ymin);
  fprintf(fpPlot, "/xmax %i def\n", xmax);
  fprintf(fpPlot, "/ymax %i def\n", ymax);
  fprintf(fpPlot, "initScaleAndFont\n");
  fprintf(fpPlot, "drawFrame\n");
  fprintf(fpPlot, "5 xticks stroke\n");
  fprintf(fpPlot, "5 yticks stroke\n");
  fprintf(fpPlot, "(%s ALPHA %.3f BETA %.3f BETA0 %.4f MAXTIME %i) title\n",
    TITLE, ALPHA, BETA, BETA0, MAXTIME);

  if (NBITER > 0) fprintf(fpPlot, "%.2f %i mt\n", TEMPS[0], COSTS[0]);
  for (i=1;i<NBITER;i++) 
    fprintf(fpPlot, "%.2f %i lt\n", TEMPS[i], COSTS[i]);

  fprintf(fpPlot, "stroke\nshowpage\n");
  fclose(fpPlot);
}
