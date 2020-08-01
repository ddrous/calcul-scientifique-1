//
// Created by drn on 10/3/19.
//

#include "neumann.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (void){

    double L = 1;

    int N = 10;     // Pas de discretidation

    neumann nmn;

    neumann_init(&nmn, L, N);

    neumann_solve(&nmn);

    neumann_display(&nmn);
}

void neumann_init (neumann *nmn, double L, int N){
    nmn->L = L;
    nmn->n = N;

    nmn->dx = L/N;

    nmn->xi = malloc(N * sizeof(double));
    nmn->u = malloc(N * sizeof(double));
    nmn->f =  malloc(N * sizeof(double));

    for (int i = 0; i < N; i++){
        nmn->xi[i] = nmn->dx * i + nmn->dx / 2;
    }

    Skyline *sky = &(nmn->sky);

    InitSkyline(sky, N);        // N est le nombre d'equations qu'on a dans notre systeme skyline a resoudre
    // InitSkyline(&(nmn->sky), N);

    // Si l'on SetSkyline alors qu'on a pas SwitchOn, alors il y a probleme  
    SwitchOn(sky, 0, 0);
    for (int i = 0; i<N-1; i++){
        SwitchOn(sky, i, i+1);
        SwitchOn(sky, i+1, i+1);        // On devrait aussi SwithOn la diagonale principale
        SwitchOn(sky, i+1, i);
    }

    AllocateSkyline(sky);

    double dx = nmn->dx;
    double val = 1/dx/dx;

    // On remplit les deux diagonales secondaires et le reste de la diagonale principale
    for(int i=0; i<N-1; i++){
        SetSkyline(sky, i, i+1, -val);
        SetSkyline(sky, i, i, 2*val + 1);
        SetSkyline(sky, i+1, i, -val);
    }
    
    // On remplit les premiers et derniers elements de la diagonales principale
    SetSkyline(sky, 0, 0, val + 1);
    SetSkyline(sky, N-1, N-1, val + 1);

    // On remplit le reste de la diagonales principale
    // for (int i = 1; i < N-1; i++){
    //     SetSkyline(sky, i, i, 2*val + 1);
    // }
    
    // Ceci display la matrice avec des etoiles a moins qu'on definisse _FULL dans la declaration de DisplaySkyline
    DisplaySkyline(sky);

    // Juste pour bien voir
    // for(int i=0; i<N-1; i++){
    //     printf("skyline (%d, %d) = %lf", i, i, GetSkyline(sky, i, i));
    //     printf("\n");
    // }

    FactoLU(sky);

    // DisplaySkyline(sky);

    // Remplissage du vecteur solution f
    for (int i = 0; i < N; i++){
        // nmn->f[i] = sin(i);
        nmn->f[i] = (double)rand()/RAND_MAX;    // La courbe de u a tendance a etre applatie aux bords peu importe le choix de f
                                                // Intuitivement, les tengeantes a u sont horizontales aux bord du domaine, d'ou l'horizontalite "presque partout"  
    }
    
}

void neumann_solve(neumann *nmn){
    SolveSkyline(&(nmn->sky), nmn->f, nmn->u);

    // Afficher la matrice M telle que MU = F
    DisplaySkyline(&(nmn->sky));

    printf("\n");
    // Afficher le deuxieme membre
    for (int i = 0; i < nmn->n; i++){
        printf("f[%d] = %f", i, nmn->f[i]);
        printf("\n");
    }

    printf("\n");
    // Afficher la solution
    for (int i = 0; i < (*nmn).n; i++){     // A noter que la notation (*nmn).n correspondant a nmn->n marche aussi tout comme en C++
        printf("u[%d] = %f", i, nmn->u[i]);
        printf("\n");
    }
    
}

void neumann_display(neumann *nmn){
    FILE *plotfile;

    plotfile = fopen("plot.dat", "w");
    for (int i = 0; i < nmn->n; i++){
        fprintf(plotfile, "%f %f %f\n", nmn->xi[i], nmn->u[i], nmn->f[i]);
    }
    
    fclose(plotfile);

    system("gnuplot plotcom");      // plotcom contient les commandes a effectuer pour afficher les graphes issues de plot.dat

}