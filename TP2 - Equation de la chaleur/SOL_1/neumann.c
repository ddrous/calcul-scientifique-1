#include "neumann.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#ifndef M_PI
  #define M_PI 3.14159265358979323846264338327950288
#endif

int main(void){


  double L = 1;

  int n = 100000;

  neumann nmn;

  neumann_init(&nmn, L, n);

  //neumann_solve(&nmn);
  double tmax = 0.001;
  serie_solve(&nmn, 1000, tmax); 
  
  neumann_display(&nmn);

}


double u0(double x){

  double u = 0;
  if (x <= 5./8 && x >= 3./8) u = 1;
  return u;

}
double f(double x){

  double y = M_PI * x;
  double v = sin(y) * sin(y) - 2 * M_PI * M_PI *
    (1 - 2 * sin(y) * sin(y));
  return v;
}

double uexact(double x){
  double y = M_PI * x;
  double v = sin(y) * sin(y);
  return v;
}

void serie_solve(neumann *nmn, int nmax,
		 double tmax){

  nmn->tmax = tmax;

  for(int i = 0; i < nmn->n; i++){
    nmn->uex[i] = 0.25;
  }
  for(int k = 1; k <= nmax;k++){
    double c = 2. / k / M_PI *
      (sin(k * M_PI * 5./8) -
       sin(k * M_PI * 3./8));
    double v = exp(-k * k * M_PI * M_PI * tmax);
    for(int i = 0; i < nmn->n; i++){
      nmn->uex[i] += c *
	 v * cos(k * M_PI * nmn->xi[i]);
    }
  }
			         
}

void neumann_init(neumann *nmn,
		  double L,
		  int n){

  nmn->L = L;
  nmn->n = n;

  nmn->dx = L/n;

  nmn->xi = malloc(n * sizeof(double));
  nmn->u = malloc(n * sizeof(double));
  nmn->uex = malloc(n * sizeof(double));
  nmn->f = malloc(n * sizeof(double));

  for(int i = 0; i < n; i++){
    nmn->xi[i] = nmn->dx * i + nmn->dx / 2;
    nmn->u[i] = u0(nmn->xi[i]);
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
    //printf("i=%d v=%f\n",i,v);
    SetSkyline(sky, i, i, v);	
  }

  for(int i = 0; i < n - 1; i++){
    v = -1 / dx /dx;
    SetSkyline(sky, i, i + 1, v);	
    SetSkyline(sky, i + 1, i, v);	
  }

  // DisplaySkyline(sky);

  
  FactoLU(sky);

  //for(int i = 0; i < n; i++)
  //  nmn->f[i] = f(nmn->xi[i]);


}


void df_solve(neumann *nmn, int pmax, double tmax){



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
    //double ue = uexact(nmn->xi[i]);
    fprintf(plotfile, "%f %f %f\n", nmn->xi[i], nmn->u[i], nmn->uex[i]);
  }

  fclose(plotfile);

  system("gnuplot plotcom");


}
