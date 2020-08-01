#include "neumann.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
int main(void){


  double L = 1;

  int n = 10;

  neumann nmn;

  neumann_init(&nmn, L, n);

  // neumann_solve(&nmn);
  double tmax = 0.001;
  df_solve(&nmn, 10, 1, tmax); 

  neumann_display(&nmn);

}

void neumann_init(neumann *nmn,
		  double L,
		  int n){

  nmn->L = L;
  nmn->n = n;

  nmn->dx = L/n;

  nmn->xi = malloc(n * sizeof(double));
  nmn->u = malloc(n * sizeof(double));
  nmn->f = malloc(n * sizeof(double));

  for(int i = 0; i < n; i++){
    nmn->xi[i] = nmn->dx * i + nmn->dx / 2;
  }

  Skyline *sky = &(nmn->sky);
  
  InitSkyline(sky, n);

  for(int i = 0; i < n - 1; i++){
    SwitchOn(sky, i, i + 1);
    SwitchOn(sky, i + 1, i);
  }

  AllocateSkyline(sky);

  double dx = nmn->dx;

  double v = 1 / dx / dx + 1;
  SetSkyline(sky, 0, 0, v);	
  SetSkyline(sky, n - 1, n - 1, v);	
  for(int i = 1; i < n - 1; i++){
    v = 2 / dx / dx + 1;
    // printf("i=%d v=%f\n",i,v);
    SetSkyline(sky, i, i, v);	
  }

  for(int i = 0; i < n - 1; i++){
    v = -1 / dx /dx;
    SetSkyline(sky, i, i + 1, v);	
    SetSkyline(sky, i + 1, i, v);	
  }

  // DisplaySkyline(sky);

  
  FactoLU(sky);

  for(int i = 0; i < n; i++)
    nmn->f[i] = 10;


}

void df_solve(neumann *nmn, int pmax, double theta, double tmax){
  Skyline M;
  Skyline N;
  InitSkyline(&M, nmn->n);
  InitSkyline(&N, nmn->n);
  
  for(int i = 0; i < nmn->n - 1; i++){
    SwitchOn(&M, i, i + 1);
    SwitchOn(&M, i + 1, i);
    SwitchOn(&N, i, i + 1);
    SwitchOn(&N, i + 1, i);
  }

  AllocateCopySkyline(&M);
  AllocateCopySkyline(&N);

  // remplir M et N .. puis faire la facto LU de chacune
  double dx = nmn->dx;
  double dt = tmax / pmax; // j'ai besoin de dt ici
  int n = nmn->n;
  
  double val = theta * dt / dx / dx + 1;
  SetSkyline(&M, 0, 0, val);	
  SetSkyline(&M, n - 1, n - 1, val);	
  for(int i = 1; i < n - 1; i++){
    val = 2 * theta * dt/ dx / dx + 1;
    // printf("i=%d val=%f\n",i,val);
    SetSkyline(&M, i, i, val);	
  }

  for(int i = 0; i < n - 1; i++){
    val = -dt / dx /dx;
    SetSkyline(&M, i, i + 1, val);	
    SetSkyline(&M, i + 1, i, val);	
  }

  DisplaySkyline(&M);

  // Remplissage de N
  val = (theta-1) * dt / dx / dx + 1;
  SetSkyline(&N, 0, 0, val);	
  SetSkyline(&N, n - 1, n - 1, val);	
  for(int i = 1; i < n - 1; i++){
    val = 2 * (theta-1) * dt/ dx / dx + 1;
    // printf("i=%d val=%f\n",i,val);
    SetSkyline(&N, i, i, val);	
  }

  for(int i = 0; i < n - 1; i++){
    val = -dt / dx /dx;
    SetSkyline(&M, i, i + 1, val);	
    SetSkyline(&M, i + 1, i, val);	
  }

  DisplaySkyline(&N);

  // FactoLU(&M);
  double v[nmn->n];

  double *u =  nmn->u;

  // double dt = tmax / pmax;
  double t =0;

  while (t<tmax){ // Voir cours 
    // v = Nu
    MatVectSkyline(&N, u, v);

    // solve Mu = v
    SolveSkyline(&M, v, u);
    t += dt;
  }

}

void neumann_solve(neumann *nmn){

  SolveSkyline(&(nmn->sky),
		   nmn->f,
		   nmn->u);

  for(int i = 0; i < nmn->n; i++){
    printf("i=%d u=%f\n",i,nmn->f[i]);
  }

}



void neumann_display(neumann *nmn){

  FILE *plotfile;

  plotfile = fopen("plot.dat", "w");

  for(int i = 0; i < nmn->n; i++){
    fprintf(plotfile, "%f %f %f\n", nmn->xi[i], nmn->u[i], nmn->u[i]);
  }

  fclose(plotfile);

  system("gnuplot plotcom");


}
