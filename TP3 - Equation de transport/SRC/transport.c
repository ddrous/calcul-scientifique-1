//
// Created by drn on 11/7/19.
// Resolution de l'equation de diffusion par la methode des differences finis, schema explicite!
//

#include "transport.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (void){

    double L = 1;
    double T = 0.5;

    int I = 100;    // Pas de discretidation en l'espace

    double cfl = 0.01;
    double c = 1;

    neumann nmn;

    // transport_init(&nmn, I, L, T, c,  cfl);
    transport_init_centre(&nmn, I, L, T, c,  cfl);
    // transport_init_lax(&nmn, I, L, T, c,  cfl);

    // transport_solve(&nmn);
    transport_solve_centre(&nmn);
    // transport_solve_lax(&nmn);

    transport_sol_exact(&nmn);

    double norme_cal = norme_1(nmn.U, nmn.dx, nmn.I);
    double norme_ext = norme_1(nmn.U_extact, nmn.dx, nmn.I);
    double erreur_rel = 100*fabs(1-norme_cal/norme_ext);
    printf("Comparaison en norme 1\n");
    printf("Solution calculee: %f\n", norme_cal);
    printf("Solution exacte: %f\n", norme_ext);
    printf("Erreur relative: %.2f\n", erreur_rel);

    norme_cal = norme_2(nmn.U, nmn.dx, nmn.I);
    norme_ext = norme_2(nmn.U_extact, nmn.dx, nmn.I);
    erreur_rel = 100*fabs(1-norme_cal/norme_ext);
    printf("\nComparaison en norme 2\n");
    printf("Solution calculee: %f\n", norme_cal);
    printf("Solution exacte: %f\n", norme_ext);
    printf("Erreur relative: %.2f\n", erreur_rel);

    norme_cal = norme_inf(nmn.U, nmn.dx, nmn.I);
    norme_ext = norme_inf(nmn.U_extact, nmn.dx, nmn.I);
    erreur_rel = 100*fabs(1-norme_cal/norme_ext);
    printf("\nComparaison en norme infinie\n");
    printf("Solution calculee: %f\n", norme_cal);
    printf("Solution exacte: %f\n", norme_ext);
    printf("Erreur relative: %.2f\n", erreur_rel);

    transport_display(&nmn);
}

double norme_1(double *U, double dx, int I){
    double norme = 0;
    for (int i = 0; i < I; i++){
        norme += dx * U[i]; 
    } 
    return norme;
}

double norme_2(double *U, double dx, int I){
    double norme = 0;
    for (int i = 0; i < I; i++){
        norme += dx * U[i]*U[i]; 
    } 
    return sqrt(norme);
}

double norme_inf(double *U, double dx, int I){
    double norme = 0;
    for (int i = 0; i < I; i++){
        if (fabs(U[i]) > norme)
            norme = fabs(U[i]); 
    } 
    return norme;
}

void transport_init(neumann *nmn, int I, double L, double T, double c, double cfl){
    nmn->I = I;
    nmn->L = L;
    nmn->dx = L/I;
    nmn->c = c;

    // nmn->dt = T/N;
    nmn->dt = cfl * nmn->dx / c ;
    nmn->T = T;
    nmn->N = T / nmn->dt;
    int N = nmn->N;

    // // Pour tester...
    // printf("\n");
    // printf("N = %d, I = %d", N, I);    
    // printf("\n");

    nmn->x_i = malloc(I * sizeof(double));
    nmn->t_n = malloc(N * sizeof(double));
    nmn->U = malloc(I * sizeof(double));
    nmn->U_0 =  malloc(I * sizeof(double));
    nmn->U_milieu =  malloc(I * sizeof(double));
    nmn->G = malloc(I * sizeof(double));

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
        SwitchOn(sky, i+1, i+1);        // On devrait aussi SwithOn la diagonale principale
        SwitchOn(sky, i+1, i);
    }

    AllocateSkyline(sky);

    double dx = nmn->dx;
    double dt = nmn->dt;
    double val = nmn->c * dt / dx;      // Valeur a utiliser pour remplir le tableau

    // On remplit les deux diagonales secondaires et le reste de la diagonale principale
    for(int i=0; i<I-1; i++){
        SetSkyline(sky, i, i, -val + 1);
        SetSkyline(sky, i + 1, i, val);
    }
    
    // On remplit les premiers et derniers elements de la diagonales principale
    SetSkyline(sky, 0, 0, -val + 1);
    SetSkyline(sky, I-1, I-1, -val + 1);
    
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
            nmn->U_0[i] = 0;    
    }
    
    
}

void transport_init_centre(neumann *nmn, int I, double L, double T, double c, double cfl){
    nmn->I = I;
    nmn->L = L;
    nmn->dx = L/I;
    nmn->c = c;

    // nmn->dt = T/N;
    nmn->dt = cfl * nmn->dx / c ;
    nmn->T = T;
    nmn->N = T / nmn->dt;
    int N = nmn->N;

    nmn->x_i = malloc(I * sizeof(double));
    nmn->t_n = malloc(N * sizeof(double));
    nmn->U = malloc(I * sizeof(double));
    nmn->U_0 =  malloc(I * sizeof(double));
    nmn->U_milieu =  malloc(I * sizeof(double));
    nmn->G = malloc(I * sizeof(double));

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
        SwitchOn(sky, i, i+1);        // On devrait aussi SwithOn la diagonale principale
        SwitchOn(sky, i+1, i);
    }

    AllocateSkyline(sky);

    double dx = nmn->dx;
    double dt = nmn->dt;
    double val = nmn->c * dt / (2*dx);      // Valeur a utiliser pour remplir le tableau

    // On remplit les deux diagonales secondaires et le reste de la diagonale principale
    for(int i=0; i<I-1; i++){
        SetSkyline(sky, i, i+1, -val);
        SetSkyline(sky, i, i, 1);
        SetSkyline(sky, i + 1, i, val);
    }
    
    // On remplit les premiers et derniers elements de la diagonales principale
    SetSkyline(sky, I-1, I-1, 1);

    for (int i = 0; i < I; i++){
            nmn->U_0[i] = 0;    
    }
      
}

void transport_init_lax(neumann *nmn, int I, double L, double T, double c, double cfl){
    nmn->I = I;
    nmn->L = L;
    nmn->dx = L/I;
    nmn->c = c;

    // nmn->dt = T/N;
    nmn->dt = cfl * nmn->dx / fabs(c) ;
    nmn->T = T;
    nmn->N = T / nmn->dt;
    int N = nmn->N;

    nmn->x_i = malloc(I * sizeof(double));
    nmn->t_n = malloc(N * sizeof(double));
    nmn->U = malloc(I * sizeof(double));
    nmn->U_0 =  malloc(I * sizeof(double));
    nmn->U_milieu =  malloc(I * sizeof(double));
    nmn->G = malloc(I * sizeof(double));

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
        SwitchOn(sky, i, i+1);        // On devrait aussi SwithOn la diagonale principale
        SwitchOn(sky, i+1, i);
    }

    AllocateSkyline(sky);

    double dx = nmn->dx;
    double dt = nmn->dt;
    double alpha = 1 - c*c * dt*dt / (dx*dx);      // Valeur a utiliser pour remplir le tableau
    double beta = (-c*dt/(2*dx)) + c*c * dt*dt / (2*dx*dx);      // Valeur a utiliser pour remplir le tableau
    double gamma = (c*dt/(2*dx)) + c*c * dt*dt / (2*dx*dx);      // Valeur a utiliser pour remplir le tableau

    // On remplit les deux diagonales secondaires et le reste de la diagonale principale
    for(int i=0; i<I-1; i++){
        SetSkyline(sky, i, i+1, beta);
        SetSkyline(sky, i, i, alpha);
        SetSkyline(sky, i+1, i, gamma);
    }
    
    // On remplit les premiers et derniers elements de la diagonales principale
    SetSkyline(sky, I-1, I-1, alpha);

    for (int i = 0; i < I; i++){
            nmn->U_0[i] = 0;    
    }
       
}

void transport_solve_lax(neumann *nmn){

    // Tableau temporaire pour stocker le resultat de la multiplication
    double * temp = malloc((nmn->I) * sizeof(double));

    // Tableau temporaire pour stocker le probleme    
    double dt = nmn->dt;
    double t = 0;
    double gamma = (nmn->c*dt/(2*nmn->dx)) + nmn->c*nmn->c * dt*dt / (2*nmn->dx*nmn->dx);      // Valeur a utiliser pour remplir le tableau

    // Tous les termes de G sont toujours nulls sauf le premier
    for (int i = 1; i < nmn->I; i++)
        nmn->G[i] = 0;

    for (int i = 0; i < nmn->I; i++)
        nmn->U[i] = nmn->U_0[i];

    while(t < nmn->T){

        // On remplit la matrice G, fonction du temps
        nmn->G[0] = exp(-t);

        MatVectSkyline(&(nmn->sky), nmn->U, temp);
        
        for (int i = 0; i < nmn->I; i++)
            nmn->U[i] = temp[i] + gamma*nmn->G[i]; 

        // Calucl de U_milieu
        if ((nmn->T/2 + dt )> t && t < (nmn->T/2 + dt)){
            for (int i = 0; i < nmn->I; i++)
                nmn->U_milieu[i] = nmn->U[i] + gamma*nmn->G[i];
        }

        t += dt;
    }
    
}

void transport_solve_centre(neumann *nmn){

    // Tableau temporaire pour stocker le resultat de la multiplication
    double * temp = malloc((nmn->I) * sizeof(double));

    // Tableau temporaire pour stocker le probleme    
    double dt = nmn->dt;
    double t = 0;
    double val = nmn->c * nmn->dt / (2*nmn->dx);

    // Tous les termes de G sont toujours nulls sauf le premier
    for (int i = 1; i < nmn->I; i++)
        nmn->G[i] = 0;

    for (int i = 0; i < nmn->I; i++)
        nmn->U[i] = nmn->U_0[i];

    while(t < nmn->T){

        // On remplit la matrice G, fonction du temps
        nmn->G[0] = exp(-t);

        MatVectSkyline(&(nmn->sky), nmn->U, temp);
        
        for (int i = 0; i < nmn->I; i++)
            nmn->U[i] = temp[i] + val*nmn->G[i]; 

        // Calucl de U_milieu
        if ((nmn->T/2 + dt )> t && t < (nmn->T/2 + dt)){
            for (int i = 0; i < nmn->I; i++)
                nmn->U_milieu[i] = nmn->U[i] + val*nmn->G[i];
        }

        t += dt;
    }
    
}

void transport_solve(neumann *nmn){

    // Tableau temporaire pour stocker le resultat de la multiplication
    double * temp = malloc((nmn->I) * sizeof(double));

    // Tableau temporaire pour stocker le probleme    
    double dt = nmn->dt;
    double t = 0;
    double val = nmn->c * nmn->dt / nmn->dx;

    // Tous les termes de G sont toujours nulls sauf le premier
    for (int i = 1; i < nmn->I; i++)
        nmn->G[i] = 0;

    for (int i = 0; i < nmn->I; i++)
        nmn->U[i] = nmn->U_0[i];

    while(t < nmn->T){

        // On remplit la matrice G, fonction du temps
        nmn->G[0] = exp(-t);

        MatVectSkyline(&(nmn->sky), nmn->U, temp);
        
        for (int i = 0; i < nmn->I; i++)
            nmn->U[i] = temp[i] + val*nmn->G[i]; 

        // Calucl de U_milieu
        if ((nmn->T/2 + dt )> t && t < (nmn->T/2 + dt)){
            for (int i = 0; i < nmn->I; i++)
                nmn->U_milieu[i] = nmn->U[i] + val*nmn->G[i];
        }

        t += dt;
    }
    

    // Afficher la matrice M telle que U_n+1 = M * U_n 
    // DisplaySkyline(&(nmn->sky));

    // printf("\n");
    // // Afficher la solution initiale
    // for (int i = 0; i < nmn->I; i++){
    //     printf("U_0[%-4d] = %f", i, nmn->U_0[i]);
    //     printf("\n");
    // }

    printf("\n");
    // // Afficher la solution
    // for (int i = 0; i < nmn->I; i++){     // A noter que la notation (*nmn).n correspondant a nmn->n marche aussi tout comme en C++
    //     printf("U_N[%-4d] = %f", i, nmn->U[i]);
    //     printf("\n");
    // }
    
}

void transport_sol_exact(neumann *nmn){
    nmn->U_extact = malloc((nmn->I)*sizeof(double));
    
    double val;
    for (int i = 0; i < nmn->I; i++){
        if(nmn->x_i[i] < nmn->c * nmn->T){
            val = -nmn->T + nmn->x_i[i]/nmn->c;
            nmn->U_extact[i] = exp(val);
        } else
            nmn->U_extact[i] = 0;
    }

    
}

void transport_display(neumann *nmn){
    FILE *plotfile;

    plotfile = fopen("plot.dat", "w");
    for (int i = 0; i < nmn->I; i++){
        fprintf(plotfile, "%f %f %f %f\n", nmn->x_i[i], nmn->U[i], nmn->U_extact[i], nmn->U_milieu[i]);
    }
    
    fclose(plotfile);

    system("gnuplot plotcom");      // plotcom contient les commandes a effectuer pour afficher les graphes issues de plot.dat

}
