#include "neumann.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
int main(void){


  double L = 1;

  int n = 10;

  neumann nmn;

  neumann_init(&nmn, L, n);

  neumann_solve(&nmn);

  neumann_display(&nmn);

}

void neumann_init(neumann *nmn, double L, int n){

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
    printf("i=%d v=%f\n",i,v);
    SetSkyline(sky, i, i, v);	
  }

  for(int i = 0; i < n - 1; i++){
    v = -1 / dx /dx;
    SetSkyline(sky, i, i + 1, v);	
    SetSkyline(sky, i + 1, i, v);	
  }

  DisplaySkyline(sky);

  
  FactoLU(sky);

  DisplaySkyline(sky);

  for(int i = 0; i < n; i++)
    nmn->f[i] = 10;


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
