//
// Created by drn on 11/7/19.
// Resolution du l'equation de transport par la methode des differences finies, theta schema
//

#include "neumann.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (void){

    double L = 1;
    double T = 0.001;

    int I = 75;    // Pas de discretidation en l'espace
    int N = 50;    // Pas de discretidation en temps

    double theta = 0.0;
    double cfl = 1.1;   // cfl < 1 pour verifier la condition CFL a la definition de dt

    neumann nmn;

    // neumann_df_explicite_init(&nmn, I, L, T);
    neumann_df_theta_init(&nmn, I, L, N, T, theta, cfl);

    // neumann_df_explicite_solve(&nmn);
    neumann_df_theta_solve(&nmn, theta);

    neumann_display(&nmn);
}

void neumann_df_explicite_init (neumann *nmn, int I, double L, double T){
    nmn->I = I;
    nmn->L = L;
    nmn->dx = L/I;

    // nmn->dt = T/nmn->N;
    nmn->dt = 0.9 * nmn->dx * nmn->dx / 2;
    nmn->T = T;
    nmn->N = T / nmn->dt;
    int N = nmn->N;

    nmn->x_i = malloc(I * sizeof(double));
    nmn->t_n = malloc(N * sizeof(double));
    nmn->U = malloc(I * sizeof(double));
    nmn->U_0 =  malloc(I * sizeof(double));

    // On remplit les differents points de discretisation en espace
    for (int i = 0; i < I; i++){
        nmn->x_i[i] = nmn->dx * i + nmn->dx / 2;
    }

    // On remplit les differents points de discretisation en temps
    for (int n = 0; n < N; n++){
        nmn->t_n[n] = nmn->dt * n;
    }

    Skyline *sky = &(nmn->sky);

    InitSkyline(sky, I);        // N est le nombre d'equations qu'on a dans notre systeme skyline a resoudre

    // Si l'on SetSkyline alors qu'on a pas SwitchOn, alors il y a probleme  
    SwitchOn(sky, 0, 0);
    for (int i = 0; i<I-1; i++){
        SwitchOn(sky, i, i+1);
        SwitchOn(sky, i+1, i+1);        // On devrait aussi SwithOn la diagonale principale
        SwitchOn(sky, i+1, i);
    }

    AllocateSkyline(sky);

    double dx = nmn->dx;
    double dt = nmn->dt;
    double val = -dt/dx/dx;      // Valeur a utiliser pour remplir le tableau

    // On remplit les deux diagonales secondaires et le reste de la diagonale principale
    for(int i=0; i<I-1; i++){
        SetSkyline(sky, i, i+1, -val);
        SetSkyline(sky, i, i, 2*val + 1);
        SetSkyline(sky, i+1, i, -val);
    }
    
    // On remplit les premiers et derniers elements de la diagonales principale
    SetSkyline(sky, 0, 0, val + 1);
    SetSkyline(sky, I-1, I-1, val + 1);
    
    // Ceci display la matrice avec des etoiles a moins qu'on definisse _FULL dans la declaration de DisplaySkyline
    // DisplaySkyline(sky);

    // Juste pour bien voir
    // for(int i=0; i<N-1; i++){
    //     printf("skyline (%d, %d) = %lf", i, i, GetSkyline(sky, i, i));
    //     printf("\n");
    // }

    // FactoLU(sky);

    // DisplaySkyline(sky);

    // Remplissage du vecteur solution initiale U_0
    for (int i = 0; i < I; i++){
        if(i > 3*(I-1)/8 && i < 5*(I-1)/8)
            nmn->U_0[i] = 1;
        else 
            nmn->U_0[i] = 0;    
    }
    
}

void neumann_df_theta_init (neumann *nmn, int I, double L, int N, double T, double theta, double cfl){
    nmn->I = I;
    nmn->L = L;
    nmn->dx = L/I;

    nmn->T = T;
    nmn->N = N;
    nmn->dt = T/N;
    // Condition CFL(Courant–Friedrichs–Lewy) pour theta < 0.5 (inclus le schema explicite theta = 0)
    // if (theta < 0.5)
        nmn->dt = cfl * nmn->dx * nmn->dx / (2 * (1-2*theta));  

    nmn->x_i = malloc(I * sizeof(double));
    nmn->t_n = malloc(N * sizeof(double));
    nmn->U = malloc(I * sizeof(double));
    nmn->U_0 =  malloc(I * sizeof(double));

    // On remplit les differents points de discretisation en espace
    for (int i = 0; i < I; i++){
        nmn->x_i[i] = nmn->dx * i + nmn->dx / 2;
    }

    // On remplit les differents points de discretisation en temps
    for (int n = 0; n < N; n++){
        nmn->t_n[n] = nmn->dt * n;
    }

    // Remplissage du vecteur solution initiale U_0
    for (int i = 0; i < I; i++){
        if(i > 3*(I-1)/8 && i < 5*(I-1)/8)
            nmn->U_0[i] = 1;
        else 
            nmn->U_0[i] = 0;    
    }
    
}

void neumann_df_explicite_solve(neumann *nmn){

    // Tableau temporaire pour stocker le resultat de la multiplication
    double * temp = malloc((nmn->I) * sizeof(double));

    MatVectSkyline(&(nmn->sky), nmn->U_0, nmn->U);

    for (int n = 0; n < nmn->N; n++){
        MatVectSkyline(&(nmn->sky), nmn->U, temp);
        
        for (int i = 0; i < nmn->I; i++)
            nmn->U[i] = temp[i];        

    }
    

    // Afficher la matrice M telle que U_n+1 = M * U_n 
    // DisplaySkyline(&(nmn->sky));

    // printf("\n");
    // // Afficher la solution initiale
    // for (int i = 0; i < nmn->I; i++){
    //     printf("U_0[%-4d] = %f", i, nmn->U_0[i]);
    //     printf("\n");
    // }

    // printf("\n");
    // // Afficher la solution
    // for (int i = 0; i < nmn->I; i++){     // A noter que la notation (*nmn).n correspondant a nmn->n marche aussi tout comme en C++
    //     printf("U_N[%-4d] = %f", i, nmn->U[i]);
    //     printf("\n");
    // }
    
}

void neumann_df_theta_solve(neumann *nmn, double theta){
    // Mauvaise pratique de definir M ici. On aurait du definir M et N dans neumann
    Skyline M;      
    Skyline N;

    // Taille des vecteurs, des matrices, etc...
    int I = nmn->I;

    InitSkyline(&M, I);   
    InitSkyline(&N, I);

    // Identification des cefficients non nules de chaque matrice
    SwitchOn(&M, 0, 0);
    SwitchOn(&N, 0, 0);
    for (int i = 0; i<I-1; i++){
        SwitchOn(&M, i, i+1);
        SwitchOn(&M, i+1, i+1);
        SwitchOn(&M, i+1, i);

        SwitchOn(&N, i, i+1);
        SwitchOn(&N, i+1, i+1);
        SwitchOn(&N, i+1, i);
    }

    AllocateSkyline(&M);
    AllocateSkyline(&N);

    double dx = nmn->dx;
    double dt = nmn->dt;


    // Remplissage de M
    double val = theta*dt/dx/dx;  
    for(int i=0; i<I-1; i++){
        SetSkyline(&M, i, i+1, -val);
        SetSkyline(&M, i, i, 2*val + 1);
        SetSkyline(&M, i+1, i, -val);
    }
    // Remplissage du premier et dernier element de la diagonle de M
    SetSkyline(&M, 0, 0, val + 1);
    SetSkyline(&M, I-1, I-1, val + 1);

    // Remplissage de N
    val = -(1-theta)*dt/dx/dx;
    for(int i=0; i<I-1; i++){
        SetSkyline(&N, i, i+1, -val);
        SetSkyline(&N, i, i, 2*val + 1);        
        SetSkyline(&N, i+1, i, -val);
    }
    SetSkyline(&N, 0, 0, val + 1);
    SetSkyline(&N, I-1, I-1, val + 1);
    
    
    // Ceci display la matrice avec des etoiles a moins qu'on definisse _FULL dans la declaration de DisplaySkyline
    // DisplaySkyline(&M);
    // DisplaySkyline(&N);

    // Juste pour bien voir
    // for(int i=0; i<N-1; i++){
    //     printf("skyline (%d, %d) = %lf", i, i, GetSkyline(&M, i, i));
    //     printf("\n");
    // }

    FactoLU(&M);

    // A l'instant initial, U = U_0
    for (int i = 0; i < I; i++){
        nmn->U[i] = nmn->U_0[i];
    }
    
    double * V = malloc(I * sizeof(double));
    double t = 0;
    while (t < nmn->T){
        MatVectSkyline(&N, nmn->U, V);
        SolveSkyline(&M, V, nmn->U);

        t += dt; 
    }
}

void neumann_display(neumann *nmn){
    FILE *plotfile;

    plotfile = fopen("plot.dat", "w");
    for (int i = 0; i < nmn->I; i++){
        fprintf(plotfile, "%f %f %f\n", nmn->x_i[i], nmn->U_0[i], nmn->U[i]);
    }
    
    fclose(plotfile);

    system("gnuplot plotcom");      // plotcom contient les commandes a effectuer pour afficher les graphes issues de plot.dat

}