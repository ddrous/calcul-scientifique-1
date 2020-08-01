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
    double  *x_i;
    double *t_n;

    // Solutions initiale, au milieu de la propagation, et solution a calculer
    double *U_0;
    double *U_milieu;
    double *U;
    double *U_extact;
    
    // Vecteur tel que U = AU + G
    double *G;

    // Vitesse de l'equation
    double c;

    // Matrice creuse du systeme (U_n+1 = sky * U_n)
    Skyline sky;

} neumann;

//Juste les signatures des fonctions a definir dans neumann.c
void transport_init(neumann *nmn, int I, double L, double T, double c, double cfl); //Afin d'eviter les copies
void transport_solve(neumann *nmn);
void transport_display(neumann *nmn);
void transport_sol_exact(neumann *nmn);

double norme_1(double *U, double dx, int I);
double norme_2(double *U, double dx, int I);
double norme_inf(double *U, double dx, int I);

void transport_init_centre(neumann *nmn, int I, double L, double T, double c, double cfl); //Afin d'eviter les copies
void transport_solve_centre(neumann *nmn);

void transport_init_lax(neumann *nmn, int I, double L, double T, double c, double cfl); //Afin d'eviter les copies
void transport_solve_lax(neumann *nmn);

#endif //SKYLINE_NEUMANN_H
