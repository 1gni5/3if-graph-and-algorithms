#include <math.h>
#include <time.h>
#include "sa.h"

void rndMove();

/* Représentation d'un point */
typedef struct {
    double x;
    double y;
}point_t;

/* Liste des points du graohe */
point_t *points;

/* Chemin dans le graphes */
int* permutation;
int* best;

/* Nombre de sommets */
int nb_points;

int swap_a, swap_b;

int cost() {
    double sum = 0;

    /* Parcours la permutation */
    for (int i = 0; i < nb_points; i++) {
        sum += hypot(
            points[permutation[(i + 1) % nb_points]].x - points[permutation[i]].x,
            points[permutation[(i + 1) % nb_points]].y - points[permutation[i]].y
        );
    }

    return sum;
}

void printConfig() {
  for (int i=0;i<nb_points;i++) printf("%i ", permutation[i]);
  printf("\n");
}

int scanConfig() {
    /* Une configuration a un nom */
    if (!scanConfigTitle()) return 0;

    /* On récupère le nombre de sommets */
    if (scanf("%i", &nb_points) != 1) {
        fprintf(stderr, "error scan number of edges\n");
        exit(1);
    }

    /* On alloue la mémoire pour les points */
    points = (point_t *)malloc(nb_points * sizeof(point_t));

    /* On récupère les points */
    for (int i = 0; i < nb_points; i++) {
        if (scanf("%lf %lf", &points[i].x, &points[i].y) != 2) {
            fprintf(stderr, "error scan point %i\n", i);
            exit(1);
        }
    }

    /* On alloue la mémoire pour la permutation */
    permutation = (int *)malloc(nb_points * sizeof(int));
    for (int i = 0; i < nb_points; i++) permutation[i] = i;
    //TODO: initialiser la permutation de manière aléatoire

    for (int i = 0; i < nb_points; i++) {
        rndMove();
    }

    /* On alloue la mémoire pour le meilleur chemin */
    best = (int *)malloc(nb_points * sizeof(int));

    scanf("\n"); /* ready to scan the next config line */
    allocTempsAndCosts();
    return 1;
}

void doMove(int i, int j) {
    int t = permutation[i];
    permutation[i] = permutation[j];
    permutation[j] = t;
}

void rndMove() {
    swap_a = rand() % nb_points;
    swap_b = rand() % nb_points;
    while (swap_a == swap_b) swap_b = rand() % nb_points;
    doMove(swap_a, swap_b);
}

void undoMove() {
    doMove(swap_b, swap_a);
}

void saveNewBest() {
    /* TODO: Utiliser memcpy */
    for (int i = 0; i < nb_points; i++) best[i] = permutation[i];
}

int
main() {
  srand(time(0)); scanParam(); initPlot();
  while (scanConfig()) {heat(); sa(); plot(); printConfig();}
  exit(0);
}