#ifndef _BURG_H
#define _BURG_H


typedef struct df_burg{

  int nx;
  double xmin, xmax;
  double Lx;

  double dx, dt, tmax, tnow;
  double cfl;

  double *xx;

  double *un;
  double *unp1;


} df_burg;

extern df_burg null_dft;

void init_df_burg(df_burg *dfc, int nx,
		  double xmin, double xmax);

void u_exact(double x, double t, double *u);

void plot_data(df_burg *dfc);

void compute_df(df_burg *dfc, double tmax, double cfl);


#endif
