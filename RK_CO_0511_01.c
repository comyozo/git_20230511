#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//
#define PI 3.141592653589793
#define PI2 2.0 * 3.141592653589793
#define ALPHA 1.0
#define BETA 1.0
#define VPL 0.1
#define FRICMAX 1.0
#define KKZ 1.0
//NEVER MOVE NEGATIVE DIRECTION
#define BK1DFOREL(xx0,xx1,xx2,xxp)  (  (xx1 + xx2 - 2.0*xx0) + KKZ*(xxp-xx0) )
#define FLAG00(xxx)   ( (int)( xxx*1.0E-10 + 9.99999999999999E-1 ) ) // FLAG of xxx >0 (->1) or <=0 (->0)
#define FLAGABTH(xxx,xth)   ( (int)( (fabs(xxx)-xth)*1.0E-10 + 9.99999999999999E-1 ) ) // FLAG of abs(xxx) - xth  > 0 (->1) or <=0 (->0)
#define FLAGAB00(xxx)       ( (int)(  fabs(xxx)     *1.0E-10 + 9.99999999999999E-1 ) ) // FLAG of abs(xxx)   > 0 (->1) or <=0 (->0)
#define FLAGSLIP(forcel,maxstfric,vvv)  ( 1- ( 1-FLAGABTH(forcel,maxstfric) )*( 1-FLAGAB00(vvv) ) ) // FLAG of ( ( abs(forcel) - maxstfric  > 0 ) OR (vvv > 0) ) (->1)
#define BK1DFORKF(vvv)  ( - ( FLAG00(vvv)-FLAG00(-vvv) ) / ( 1.0+ALPHA*fabs(vvv) ) ) // Kinetic Friction Force
#define ReLU(vvv)  ( 0.5*( vvv + fabs(vvv) ) )
#define FFVV(forcel,maxstfric,vvv)    ( FLAGSLIP(forcel,maxstfric,vvv)*( forcel + BETA*BK1DFORKF(vvv) ) )
//
#define NN 2
//#define KK 1.0
#define OMEGA 1.0
#define FILENAME "dat0522-01.dat"

//
void rk4(double, double, double *, double *, double *, double *);
void derivxx(double, double *, double *, double *);
void derivyy(double, double *, double *, double *);
//

int main(void)
{
  double dt, maxtime = 100.0,time;
  double xx[NN + 2], yy[NN + 2], xxout[NN + 2], yyout[NN + 2];
  FILE *fpcur;
  int imaxtime;

  fpcur = fopen(FILENAME, "w");

  /* initial condition*/
  xx[1] = -1.0;
  xx[2] = 1.0;
  yy[1] =  0.0;
  yy[2] =  0.0;
  
 /* initial condition*/
 
  dt = 0.01;
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

    printf("time = %10.3lf xx1= %10.3lf   xx2=  %10.3lf  \n", time, xx[1], xx[2]);
    fprintf(fpcur, "%10.3lf %10.3lf %10.3lf %10.3lf %10.3lf \n", time, xx[1], yy[1],xx[2], yy[2]);
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
    yytt[ii] = ReLU(yy[ii] + dthh * ddyy1[ii]);
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
    yytt[ii] = ReLU(yy[ii] + dthh * ddyy2[ii]);
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
    yytt[ii] = ReLU(yy[ii] + dt * ddyy3[ii]);
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
    yyout[ii] = ReLU(yy[ii] + dt6 * (ddyy1[ii] + 2.0 * ddyy2[ii] + 2.0 * ddyy3[ii] + ddyy4[ii]));
  }
  /* end of 4-thd step*/
}
//
void derivxx(double time, double *xx, double *yy, double *ddxx)
{
  for (int ii = 1; ii <= NN; ii++)    ddxx[ii] = ReLU(yy[ii]);
}
//
void derivyy(double time, double *xx, double *yy, double *ddyy)
{
  double forceel,xpph;

  xpph=VPL*time;

  /*for (int ii = 1; ii <= NN; ii++)    ddyy[ii] = xx[ii + 1] + xx[ii - 1] - 2.0 * xx[ii]+KK*(time-xx[ii]);*/
  for(int ii=1;ii<=NN;ii++)
  {
    forceel=BK1DFOREL(xx[ii],xx[ii-1],xx[ii+1],xpph);
    ddyy[ii]=FFVV(forceel,FRICMAX,yy[ii]);
  }
}
