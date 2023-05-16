#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//
#define PI 3.141592653589793
#define PI2 2.0 * 3.141592653589793
#define NN 4
//#define KK 1.0
#define OMEGA 1.0
#define FILENAME "dat3.dat"

//
void rk4(double, double, double *, double *, double *, double *);
void derivxx(double, double *, double *, double *);
void derivyy(double, double *, double *, double *);
//
void main()
{
  double dt, maxtime = 20.0,time;
  double xx[NN + 2], yy[NN + 2], xxout[NN + 2], yyout[NN + 2];
  FILE *fpcur;
  int imaxtime;

  fpcur = fopen(FILENAME, "w");

  /* initial condition*/
  xx[1] = 1.0;
  xx[2] = 0.0;
  xx[3] = 0.0;
  xx[4] = -1.0;
  yy[1] =  0.0;
  yy[2] =  0.0;
  yy[3] =  0.0;
  yy[4] =  0.0;
 /* initial condition*/
 
  dt = 0.01*PI2/OMEGA;
  imaxtime=maxtime/dt;

  for (int itime = 0; itime <= imaxtime; itime++)
  {
    
    time = dt * itime;

    rk4(dt, time, xx, yy, xxout, yyout);

    for (int ii = 1; ii <= NN; ii++) 
    {
      xx[ii] = xxout[ii];
      yy[ii] = yyout[ii];
    }

    printf("time = %10.3lf xx1= %10.3lf   xx2=  %10.3lf  xx3=  %10.3lf xx4=  %10.3lf \n", time, xx[1], xx[2],xx[3],xx[4]);
    fprintf(fpcur, "%10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf\n", time, xx[1], yy[1],xx[2], yy[2],xx[3], yy[3],xx[4], yy[4]);
  }

  fclose(fpcur);
}

void rk4(double dt, double time, double *xx, double *yy, double *xxout, double *yyout)
{
  double xxtt[NN + 2], yytt[NN + 2];
  double ddxx1[NN + 2], ddyy1[NN + 2], ddxx2[NN + 2], ddyy2[NN + 2], ddxx3[NN + 2], ddyy3[NN + 2], ddxx4[NN + 2], ddyy4[NN + 2];
  double dthh = dt / 2.0, dt6 = dt / 6.0;

  xx[0]      = xx[1];
  xx[NN + 1] = xx[NN];

  /* start of 1st step*/
  derivxx(time, xx, yy, ddxx1);
  derivyy(time, xx, yy, ddyy1);
  for (int ii = 1; ii <= NN; ii++)
  {
    xxtt[ii] = xx[ii] + dthh * ddxx1[ii];
    yytt[ii] = yy[ii] + dthh * ddyy1[ii];
  }
  //
  xxtt[0]      = xxtt[1];
  xxtt[NN + 1] = xxtt[NN];
  /* end of 1st step*/
  //
  /* start of 2nd step*/
  derivxx(time + dthh, xxtt, yytt, ddxx2);
  derivyy(time + dthh, xxtt, yytt, ddyy2);
  for (int ii = 1; ii <= NN; ii++)
  {
    xxtt[ii] = xx[ii] + dthh * ddxx2[ii];
    yytt[ii] = yy[ii] + dthh * ddyy2[ii];
  }
  //
  xxtt[0]      = xxtt[1];
  xxtt[NN + 1] = xxtt[NN];
  /* end of 2nd step*/
  //
  /* start of 3nd step*/
  derivxx(time + dthh, xxtt, yytt, ddxx3);
  derivyy(time + dthh, xxtt, yytt, ddyy3);
  for (int ii = 1; ii <= NN; ii++)
  {
    xxtt[ii] = xx[ii] + dt * ddxx3[ii];
    yytt[ii] = yy[ii] + dt * ddyy3[ii];
  }
  //
  xxtt[0]      = xxtt[1];
  xxtt[NN + 1] = xxtt[NN];
  /* end of 3nd step*/
  //
  /* start of 4th step*/
  derivxx(time + dt, xxtt, yytt, ddxx4);
  derivyy(time + dt, xxtt, yytt, ddyy4);
  for (int ii = 1; ii <= NN; ii++)
  {
    xxout[ii] = xx[ii] + dt6 * (ddxx1[ii] + 2.0 * ddxx2[ii] + 2.0 * ddxx3[ii] + ddxx4[ii]);
    yyout[ii] = yy[ii] + dt6 * (ddyy1[ii] + 2.0 * ddyy2[ii] + 2.0 * ddyy3[ii] + ddyy4[ii]);
  }
  /* end of 4-thd step*/
}
//
void derivxx(double time, double *xx, double *yy, double *ddxx)
{
  for (int ii = 1; ii <= NN; ii++)    ddxx[ii] = yy[ii];
}
//
void derivyy(double time, double *xx, double *yy, double *ddyy)
{
  for (int ii = 1; ii <= NN; ii++)    ddyy[ii] = xx[ii + 1] + xx[ii - 1] - 2.0 * xx[ii];
}