#ifndef _NEUMANN_H
#define _NEUMANN_H

#include "skyline.h"

typedef struct neumann{

  // nombre de points
  int n;
  double L;
  double dx;
  double tmax;

  // points du maillages
  double *xi;

  // solution et second membre
  double *u;
  double *f;

  // solution "exacte"
  double *uex;
  
  // matrice creuse
  Skyline sky;


} neumann;

void neumann_init(neumann *nmn,
		  double L,
		  int n);

void serie_solve(neumann *nmn, int nmax, double tmax);
void df_solve(neumann *nmn, int pmax, double tmax);
void neumann_solve(neumann *nmn);
void neumann_display(neumann *nmn);



#endif
