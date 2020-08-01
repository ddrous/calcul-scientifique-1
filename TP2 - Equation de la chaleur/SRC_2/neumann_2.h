#ifndef _NEUMANN_H
#define _NEUMANN_H

#include "skyline.h"

typedef struct neumann{

  // nombre de points
  int n;
  double L;
  double dx;

  // points du maillages
  double *xi;

  // solution et second membre
  double *u;
  double *f;

  // matrice creuse
  Skyline sky;


} neumann;

void neumann_init(neumann *nmn,
		  double L,
		  int n);
void neumann_solve(neumann *nmn);
void df_solve(neumann *nmn, int pmax, double theta, double tmax);
void neumann_display(neumann *nmn);



#endif
