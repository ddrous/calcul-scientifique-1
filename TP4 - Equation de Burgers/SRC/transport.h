//
// Created by drn on 10/3/19.
// Resolution du l'equation de la chaleur pa la methode des differences finies, schema explicite 
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
    int N;      // Nombre de ponts du maillage en temps
    double T;   // Longueur du maillage en temps
    double dt;

    //  Differents points du maillage en temops et en espace
    double x_min;
    double x_max;
    double  *x_i;
    double *t_n;

    // Solutions initiale, au milieu de la propagation, et solution a calculer
    double *U_0;
    double *U_milieu;
    double *U;
    double *U_extact;

    double epsilon;

} neumann;

//Juste les signatures des fonctions a definir dans neumann.c
void transport_init(neumann *nmn, int I, double L, int N, double T, double x_min, double x_max, double epsilon, double cfl); //Afin d'eviter les copies
void transport_solve(neumann *nmn, double cfl);
void transport_solve_trafic(neumann *nmn, double cfl);
void transport_solve_df(neumann *nmn);
void transport_display_trafic(neumann *nmn);
void transport_display(neumann *nmn);
void transport_sol_exact(neumann *nmn);

double f(double x);
double f_prime(double x);
double F(double x, double y, double lamda);


#endif //SKYLINE_NEUMANN_H
