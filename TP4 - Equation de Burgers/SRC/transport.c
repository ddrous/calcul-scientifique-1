//
// Created by drn on 11/7/19.
// Resolution de l'equation de diffusion par la methode des differences finis, schema explicite!
//

#include "transport.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double norme_sup(double *U, int N){
    double norme = 0;
    for (size_t i = 0; i < N; i++){
        if (fabs(U[i]) > norme)
            norme = fabs(U[i]);
    }
    return norme;
}

int main (void){

    double T = 2;

    int I = 10000;    // Pas de discretidation en l'espace
    double N = 100;     // N'influence pas le resultat

    double x_min = -1;
    // double x_min = 0;
    double x_max = 2;
    // double x_max = 1;
    double L = x_max-x_min;
    double epsilon = 0.1;
    double cfl = 0.5;

    neumann nmn;

    transport_init(&nmn, I, L, N, T, x_min, x_max, epsilon, cfl);

    transport_sol_exact(&nmn);

    // transport_solve(&nmn, cfl);

    transport_solve_df(&nmn);


    // printf("f(%f) = %f\n", 5., f(5));


    transport_display(&nmn);

    // FILE *plotfile;
    // plotfile = fopen("data.trafic.dat", "w");
    // double T = 0;
    // while (true){
    //     // neumann nmn;
    //     transport_init(&nmn, I, L, N, T, x_min, x_max, epsilon, cfl);
    //     transport_solve_trafic(&nmn, cfl);
    //     fprintf(plotfile, "%-8f %-8f\n", T, norme_sup(nmn.U, nmn.I));
    //     if (norme_sup(nmn.U, nmn.I) < 1e-10)
    //         break;
    //     T += 0.05;
    // }
    // fclose(plotfile);
    // system("gnuplot plotcom_trafic_densite");
    // printf("Le temps pour se vider est de: %f\n", T);
    
    // transport_display_trafic(&nmn);


}

double f(double x){
    return x*x/2;
    // return 130*x*(1-x/200);
}

double f_trafic(double x){
    return 130*x*(1-x/200);
}

double f_prime(double x){
    return x;
    // return 130*(1-2*x/200);
}

double f_prime_trafic(double x){
    return 130*(1-2*x/200);
}

double F(double x, double y, double lamda){
    return (f(x)+f(y)-lamda*(y-x))/2;
}

double F_trafic(double x, double y, double lamda){
    return (f_trafic(x)+f_trafic(y)-lamda*(y-x))/2;
}


void transport_init(neumann *nmn, int I, double L, int N, double T, double x_min, double x_max, double epsilon, double cfl){
    nmn->I = I;
    nmn->L = L;
    nmn->x_min = x_min;
    nmn->x_max = x_max;
    nmn->dx = (x_max-x_min)/I;
    nmn->epsilon = epsilon;

    // nmn->dt = T/N;
    double lambda = 2*nmn->epsilon / nmn->dx;
    double v = lambda / 2;      // car surp(f') = 2
    // nmn->dt = cfl*nmn->dx/v;
    nmn->dt = cfl*nmn->dx;
    // nmn->dt = T/N;
    nmn->T = T;
    nmn->N = N;

    // printf("dt direct = %f, avec lambda = %f\n", T/N, nmn->dt);
    // printf("epsilon = %f\n", epsilon);
    // printf("dx = %f\n", nmn->dx);
    // printf("lambda = %f\n", 2*nmn->epsilon / nmn->dx);
    // printf("v = %f\n", v);

    nmn->x_i = malloc(I * sizeof(double));
    nmn->t_n = malloc(N * sizeof(double));
    nmn->U = malloc((I+1) * sizeof(double));
    nmn->U_0 =  malloc(I * sizeof(double));
    nmn->U_milieu =  malloc(I * sizeof(double));
    nmn->U_extact = malloc((nmn->I)*sizeof(double));

    // On remplit les differents points de discretisation en espace
    for (int i = 0; i < I; i++){
        nmn->x_i[i] = nmn->x_min + nmn->dx * i + nmn->dx / 2;
    }

    // On remplit les differents points de discretisation en temps
    for (int n = 0; n < N; n++){
        nmn->t_n[n] = nmn->dt * n;
    }

    
    // Remplissage du vecteur solution initiale U_0
    for (int i = 0; i < I; i++){
            if (nmn->x_i[i] < 0)
                nmn->U_0[i] = 1;
            else{ if (nmn->x_i[i] <= 1)
                nmn->U_0[i] = 1-nmn->x_i[i];    
            else
                nmn->U_0[i] = 0;
            }
    }  

    // for (int i = 0; i < I; i++){
    //     printf("U[%d] = %f\n", i, nmn->U_0[i]);
    // }
      
}


void transport_solve(neumann *nmn, double cfl){

    double *temp = malloc((nmn->I+1)* sizeof(double));

    // Tableau temporaire pour stocker le probleme    
    double dt = nmn->dt;
    double t = 0;
    double lambda = 2*nmn->epsilon / nmn->dx;
    double lambda_max = 0;

    for (int i = 0; i < nmn->I; i++)
        nmn->U[i] = nmn->U_0[i];

    while(t < nmn->T){
    // for (int i = 0; i < nmn->N; i++){
    
        bool break_val = false;
        // MatVectSkyline(&(nmn->sky), nmn->U, temp);
        nmn->U[0] = 1; 
        
        for (int i = 0; i < nmn->I; i++)
            temp[i] = nmn->U[i];

        temp[nmn->I] = 0;
        // temp[nmn->I] = nmn->U[0];

        lambda_max = 0;
        for (int i = 0; i < nmn->I; i++){
            if(fabs(f_prime(temp[i])) > lambda_max)
                lambda_max = fabs(f_prime(temp[i]));
        }

        for (int i = 1; i < nmn->I; i++){

            if(isnan(nmn->U[i]) || isnan(nmn->U[i+1]) || isnan(nmn->U[i-1]))
                continue;

            if(fabs(f_prime(temp[i+1]) > fabs(f_prime(temp[i-1]))))
                lambda = fabs(f_prime(temp[i+1]));
            else
                lambda = fabs(f_prime(temp[i-1]));

            // lambda = fabs((f(temp[i+1]) - f(temp[i-1])) / (temp[i+1] - (temp[i-1])));

            nmn->dt = cfl * nmn->dx / lambda_max;        // Car CFL ==> lambda <= dx/dt
            

            // nmn->dt = nmn->T/nmn->N;
            // nmn->dx = nmn->dt * lambda / cfl;

            if(i==1 && t==1*dt){
                printf("lambda*dx = %f\n", lambda*nmn->dx);
                printf("at the end lambda = %f, temp[i+1] = %f, temp[i-1] = %f\n", lambda, temp[i+1], temp[i-1]);
            }


            // nmn->U[i] = temp[i] + nmn->dt * (f(temp[i+1])-f(temp[i-1])/(2*nmn->dx))
            //                     - lambda * nmn->dt * (temp[i+1]-2*temp[i]+temp[i-1])/ (2*nmn->dx); 


            nmn->U[i] = temp[i] - nmn->dt *(F(temp[i], temp[i+1], lambda)-F(temp[i-1], temp[i], lambda)) / nmn->dx; 

            if(isnan(nmn->U[i])){
                // break_val = true;
                // break;
            }
            
            // if (t == nmn->T - dt){
            //     printf("F(%f) = %f\n",temp[i], F(temp[i], temp[i], lambda));
            // }

        }

        if (break_val){
            printf("break at t = %f\n", t);
            break;
            // continue;
        }

        // printf("at the end lambda = %f, temp[i+1] = %f\n", lambda, temp[nmn->I-1]);

        t += dt;
    }
    
}

void transport_solve_trafic(neumann *nmn, double cfl){

    double *temp = malloc((nmn->I+1)* sizeof(double));
    // Tableau temporaire pour stocker le probleme    
    double dt = nmn->dt;
    double t = 0;
    double lambda = 2*nmn->epsilon / nmn->dx;
    double lambda_max = 0;

    for (int i = 0; i < nmn->I; i++)
        nmn->U[i] = 0;

    for (int i = 1; i < nmn->I/1; i++){
        // nmn->U[i] = 0.25*nmn->x_i[i]*(50-170) + 170;
        // nmn->U[i] = 2*nmn->x_i[i]*(50-170) + 170;
        nmn->U[i] = 200;
    }
    
    // for (int i = 0; i < nmn->I; i++)
    //     nmn->U[i] = 0;

    while(t < nmn->T){
        bool break_val = false;
        
        temp[0] = 0;
        for (int i = 1; i < nmn->I; i++)
            temp[i] = nmn->U[i];
        temp[nmn->I-1] = 0;


        lambda_max = 0;
        for (int i = 0; i < nmn->I; i++){
            if(fabs(f_prime_trafic(temp[i])) > lambda_max)
                lambda_max = fabs(f_prime_trafic(temp[i]));
        }

        for (int i = 1; i < nmn->I; i++){

            if(isnan(nmn->U[i]) || isnan(nmn->U[i+1]) || isnan(nmn->U[i-1]))
                continue;

            if(fabs(f_prime(temp[i+1]) > fabs(f_prime(temp[i-1]))))
                lambda = fabs(f_prime_trafic(temp[i+1]));
            else
                lambda = fabs(f_prime_trafic(temp[i-1]));

            nmn->dt = cfl * nmn->dx / lambda_max;        // Car CFL ==> lambda <= dx/dt

            // if(i==1 && t==1*dt){
            //     printf("lambda*dx = %f\n", lambda*nmn->dx);
            //     printf("at the end lambda = %f, temp[i+1] = %f, temp[i-1] = %f\n", lambda, temp[i+1], temp[i-1]);
            // }

            nmn->U[i] = temp[i] - nmn->dt *(F_trafic(temp[i], temp[i+1], lambda)-F_trafic(temp[i-1], temp[i], lambda)) / nmn->dx; 

            if(isnan(nmn->U[i])){
                // break_val = true;
                // break;
            }
        }
        if (break_val){
            printf("break at t = %f\n", t);
            break;
            // continue;
        }
        t += dt;
    }

};


void transport_solve_df(neumann *nmn){

    double *temp = malloc((nmn->I+1)* sizeof(double));

    // Tableau temporaire pour stocker le probleme    
    double dt = nmn->dt;
    double t = 0;
    double lambda = 2*nmn->epsilon / nmn->dx;


    for (int i = 0; i < nmn->I; i++)
        nmn->U[i] = nmn->U_0[i];

    while(t < nmn->T){

        // MatVectSkyline(&(nmn->sky), nmn->U, temp);
        nmn->U[0] = 1; 
        
        for (int i = 0; i < nmn->I; i++)
            temp[i] = nmn->U[i];

        temp[nmn->I] = 0;


        for (int i = 1; i < nmn->I; i++){
            nmn->U[i] = temp[i] - nmn->dt*temp[i-1]*(temp[i]-temp[i-1])/(nmn->dx); 

        }
        t += dt;
    }

    
}

void transport_sol_exact(neumann *nmn){
    
    double val;
    for (int i = 0; i < nmn->I; i++){
            if (nmn->x_i[i] < nmn->T){
                nmn->U_extact[i] = 1;
            }else{ 
                if (nmn->x_i[i] <= 1)
                    nmn->U_extact[i] = (1-nmn->x_i[i])/(1-nmn->T);    
                else
                    nmn->U_extact[i] = 0;
            }

            if(nmn->T > 1)
                nmn->U_extact[i] = sqrt(-1);    // la solution n'existe pas
    }

    
}

void transport_display(neumann *nmn){
    FILE *plotfile;

    plotfile = fopen("plot.dat", "w");
    for (int i = 0; i < nmn->I; i++){
        fprintf(plotfile, "%-8f %-8f %-8f %-8f\n", nmn->x_i[i], nmn->U[i], nmn->U_extact[i], nmn->U_0[i]);
    }
    
    fclose(plotfile);

    system("gnuplot plotcom");      // plotcom contient les commandes a effectuer pour afficher les graphes issues de plot.dat

}

void transport_display_trafic(neumann *nmn){
    FILE *plotfile;

    plotfile = fopen("plot.dat", "w");
    for (int i = 0; i < nmn->I; i++){
        fprintf(plotfile, "%-8f %-8f %-8f %-8f\n", nmn->x_i[i], nmn->U[i], nmn->U_extact[i], nmn->U_0[i]);
    }
    
    fclose(plotfile);

    system("gnuplot plotcom_trafic");      // plotcom contient les commandes a effectuer pour afficher les graphes issues de plot.dat
};

