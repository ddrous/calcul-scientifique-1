//
// Created by drn on 10/3/19.
//

#ifndef SKYLINE_NEUMANN_H
#define SKYLINE_NEUMANN_H

#include <stdbool.h>
#include "skyline.h"

typedef struct neumann{

    // nombre de points
    int n;
    double L;
    double dx;

    //Points du maillage
    double  *xi;

    //Solutions et second membre
    double *u;
    double *f;

    //matrice creuse
    Skyline sky;

} neumann;

//Juste les signatures des fonctions a definir dans neumann.c
void neumann_init(neumann *nmm, double L, int N); //Afin d'eviter les copies
void neumann_solve(neumann *nmn);
void neumann_display(neumann *nmn);

#endif //SKYLINE_NEUMANN_H
