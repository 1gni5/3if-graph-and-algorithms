#include "sa.h"
#include <time.h>

static int N;            /* chessboard has N*N squares                       */
static int SDIAG;        /* chessboard has SDIAG diagonals                   */
static int *B;           /* B[i] is the row of the queen on col i            */
                         /* the queens are initially on the main diagonal    */
static int MVI, MVJ;     /* a move swap the rows of queens on cols MVI & MVJ */
                         /* two queens are in conflict when on a same diag   */
static int *D1;          /* main diagonals                                   */
static int *D2;          /* anti diagonals                                   */

int /* number of pairs of queens attacking each other */
cost() {
  int k, d = 0;
  for (int i=0;i<SDIAG;i++) {D1[i]=0; D2[i]=0;}
  for (int i=0;i<N;i++) {
    k = (i-B[i])+N-1; /* all squares of the main diagonal controlled by 
                         queen i have the same k                    */
    if (D1[k]==0) D1[k]++; else d++;
    k = (i+B[i]);     /* all squares of the anti diagonal controlled by 
                         queen i have the same k                    */
    if (D2[k]==0) D2[k]++; else d++;
  }
  return d;
}

void
printConfig() {
  for (int i=0;i<N;i++) printf("%i ", B[i]);
  printf("\n");
}

int
scanConfig() {
  if (!scanConfigTitle()) return 0;
  if (scanf("%i", &N) != 1) {
    fprintf(stderr, "error scan number of tiles\n");
    exit(1);
  }
  SDIAG = 2*N-1;
  free(B);
  B = (int *)malloc(N*sizeof(int));
  for (int i=0 ; i<N ; i++) B[i] = i;
  scanf("\n"); /* ready to scan the next config line */
  free(D1); D1 = (int *)malloc(SDIAG*sizeof(int));
  free(D2); D2 = (int *)malloc(SDIAG*sizeof(int));
  allocTempsAndCosts();
  return 1;
}

void /* swap rows of queens on columns i and j */
doMove(int i, int j) { int t = B[i]; B[i] = B[j]; B[j] = t; }

void
rndMove() {
  MVI = MVJ = rand() % N;
  while (MVI == MVJ) MVJ = rand() % N;
  doMove(MVI,MVJ);
}

void undoMove() { doMove(MVI,MVJ); }

void saveNewBest() {} /* no need, the final conf of cost 0 is the best one */

int
main() {
  srand(time(0)); scanParam(); initPlot();
  while (scanConfig()) {heat(); sa(); plot(); printConfig();}
  exit(0);
}
