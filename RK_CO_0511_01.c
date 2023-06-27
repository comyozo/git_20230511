#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//
#define PI 3.141592653589793
#define PI2 2.0 * 3.141592653589793
#define RAND_MAX 0x7fffffff
#define ACCURACY00 1.0E-1
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
//#define MMIN -10
//#define MMAX 10
#define MNN 10
//
#define FLAGINITRAND  1 // initial displacement is given by random numbers, otherwise initial displacement is 0.0
#define ISEED  728761879   // seed of random number
#define WIRAND 0.5   // width of initial displacement given by random number
#define EVENTNO 100000
//
#define NN 10
//#define KK 1.0
#define OMEGA 1.0
#define FILENAME "dat0620-01.dat"
//
//
void rk4(double, double, double *, double *, double *, double *);
void derivxx(double, double *, double *, double *);
void derivyy(double, double *, double *, double *);
//

int main(void)
{
  double dt, maxtime = 50000.0, time;
  double xx[NN + 2], yy[NN + 2], xxout[NN + 2], yyout[NN + 2];
  FILE *fpcur;
  int imaxtime;
  int flagslip, no_of_earthquake;
  double xxav0;
  double  eqtimestt[EVENTNO], eqtimesend[EVENTNO], eqmagnitude[EVENTNO];

  double mmax = 0;
  double mmin = 0;

  fpcur = fopen(FILENAME, "w");

//Initialization
    if( FLAGINITRAND )
    {
        srand(ISEED);
        /*printf("RAND_MAX: %d\n", RAND_MAX);*/
        for(int ii = 1; ii <= NN; ii++)
        {
            double xrand;
            xrand= WIRAND*( rand()/(double)RAND_MAX - 0.5);
            xx[ii] = xrand;
            yy[ii] = 0.0;
        }
    }
    else
    {
        for(int ii = 1; ii <= NN; ii++)
        {
            xx[ii] = 0.0;
            yy[ii] = 0.0;
        }
    }

  /*printf("RAND_MAX: %d\n", RAND_MAX);
  for(int ii = 1; ii <= NN; ii++) 
  printf("xx[%3d]=  %10.3lf\n",ii,xx[ii]);*/ 

  /* initial condition*/
  //xx[1] = -1.0;
  //xx[2] =  1.0;
  //yy[1] =  0.0;
  //yy[2] =  0.0;
  /* initial condition*/

  flagslip = 0;
  no_of_earthquake = 0;
  xxav0 = 0.0;
  for(int ii = 1; ii <= NN; ii++) 
    xxav0 += xx[ii];
 
  dt = 0.01;
  imaxtime = maxtime / dt;

  for (int itime = 0; itime <= imaxtime; itime++)
  {
    double xxav,yyav,sumslip,yymax,magnitude;
    int itime00;
    
    time = dt * itime;

    rk4(dt, time, xx, yy, xxout, yyout);

    xxav = 0.0;
    yyav = 0.0;
    yymax = 0.0;

    for(int ii = 1; ii <= NN; ii++)
    {
      xxav += xxout[ii];
      yyav += yyout[ii];
      yymax = ReLU( (yymax - yyout[ii]) ) + yyout[ii];
    }

//judgement of the beginning of slipping
        if(flagslip == 0 && ( yymax >= ACCURACY00) )
        {
           flagslip = 1;
           //sumslip = 0.0;
           //itime00 = itime;
           xxav0 = xxav;
           eqtimestt[no_of_earthquake] = time;
           printf("A. No. = %5d \n", no_of_earthquake);
        }

//judgement of the end of slipping
        if(flagslip == 1 && ( yymax < ACCURACY00) )
        {
           flagslip = 0;
           sumslip = xxav - xxav0;
           magnitude = log(sumslip);
           eqtimesend[no_of_earthquake] = time;
           eqmagnitude[no_of_earthquake] = magnitude;

           if(magnitude > mmax)
           {
              mmax = magnitude;
           }

           if(magnitude < mmin)
           {
              mmin = magnitude;
           }

           /*printf("B. No. = %5d \n", no_of_earthquake);*/
           no_of_earthquake++;
        }

    for (int ii = 1; ii <= NN; ii++) 
    {
      xx[ii] = xxout[ii];
      yy[ii] = yyout[ii];
    }

    /*printf("time = %10.3lf xx1= %10.3lf   xx2=  %10.3lf    xx3= %10.3lf   xx4=  %10.3lf\n", time, xx[1], xx[2], xx[3], xx[4]);*/
    /*fprintf(fpcur, "%10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf\n", time, xx[1], yy[1],xx[2], yy[2], xx[3], yy[3], xx[4], yy[4]);*/
  }

  /*for(int ii = 0; ii < no_of_earthquake; ii++)  
  printf("No. = %5d tstt= %10.3lf   tend=   %10.3lf  mag = %10.3lf\n", ii, eqtimestt[ii], eqtimesend[ii],eqmagnitude[ii]);
  printf("No. = %5d \n", no_of_earthquake);

  fclose(fpcur);*/

  double neqdev = ((double)(mmax - mmin)) / MNN;
  int neq[MNN];

  for (int ii = 0; ii == MNN-1; ii++)
  {
    neq[ii] = 0;
  }

  for (int ii = 0; ii <= no_of_earthquake - 1; ii++)
  {
    neq[(int)((eqmagnitude[ii] - mmin)/neqdev)]++;
  }

  for (int ii = 0; ii <= MNN; ii++)
  {
    printf("%lf %d\n", ii*neqdev + mmin, neq[ii]);
    fprintf(fpcur, "%lf %d\n", ii*neqdev + mmin, neq[ii]);
  }

  printf("mmin = %lf \nmmax = %lf \n", mmin, mmax);
  
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
    yytt[ii] = ReLU( yy[ii] + dthh * ddyy1[ii] );
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
    yytt[ii] = ReLU( yy[ii] + dthh * ddyy2[ii] );
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
    yytt[ii] = ReLU( yy[ii] + dt * ddyy3[ii] );
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
    yyout[ii] = ReLU( yy[ii] + dt6 * (ddyy1[ii] + 2.0 * ddyy2[ii] + 2.0 * ddyy3[ii] + ddyy4[ii]) );
  }
  /* end of 4-thd step*/
}
//
void derivxx(double time, double *xx, double *yy, double *ddxx)
{
  for (int ii = 1; ii <= NN; ii++)    ddxx[ii] = ReLU( yy[ii] );
}
//
void derivyy(double time, double *xx, double *yy, double *ddyy)
{
  double forceel, xpph;

  xpph=VPL*time;

  /* for (int ii = 1; ii <= NN; ii++)    ddyy[ii] = xx[ii + 1] + xx[ii - 1] - 2.0 * xx[ii];*/
  for (int ii = 1; ii <= NN; ii++)
  {
    forceel = BK1DFOREL(xx[ii], xx[ii-1], xx[ii+1], xpph);
    ddyy[ii] = FFVV(forceel, FRICMAX, yy[ii]);
  }
}
