//
// Created by drn on 10/3/19.
// Resolution du l'equation de transport pa la methode des differences finies, theta schema
//

#ifndef SKYLINE_NEUMANN_H
#define SKYLINE_NEUMANN_H

#include <stdbool.h>
#include "skyline.h"

typedef struct neumann{

    // EN ESPACE
    int I;      // Nombre de points du maillage en espace
    double L;   // Longueur du maillage en espace   
    double dx;

    // EN TEMPS
    int N;      // Nombre de points du maillage en temps
    double T;   // Longueur du maillage en temps
    double dt;

    //  Differents points du maillage en temops et en espace
    double  *x_i;
    double *t_n;

    // Solutions initiale et solution a calculer
    double *U_0;
    double *U;

    // Matrice creuse du systeme (U_n+1 = sky * U_n)
    Skyline sky;

} neumann;

//Juste les signatures des fonctions a definir dans neumann.c
void neumann_df_explicite_init(neumann *nmn, int I, double L, double T);                 // Initialisation du systeme pour le schema explicite (Voire TP2)
void neumann_df_theta_init(neumann *nmn, int I, double L, int N, double T, double theta, double cfl);      // Initialisation du systeme pour le theta schema
void neumann_df_explicite_solve(neumann *nmn);
void neumann_df_theta_solve(neumann *nmn, double theta);                  // Resolotion du systeme pour le theta schema
void neumann_display(neumann *nmn);

#endif //SKYLINE_NEUMANN_H
