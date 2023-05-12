#include "D_Cont.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define qave_EXP  0.547
#define qave_COEF 0.386
#define qdev_EXP  0.673
#define qdev_COEF 1.327
#define L_max 2000.0
//#define A 0.0 //5.0
int Ini_Dcont=0;
float D0;

float Compute_qinf(int L1, int L2, int NORM)
{
  double L12;

  if(NORM==0){if(L1<L2){L12=L1;}else{L12=L2;}}
  else if(NORM==1){if(L1>L2){L12=L1;}else{L12=L2;}}
  else{L12=sqrt((float)(L1*L2));}

  double qave=qave_COEF*pow(L12, -qave_EXP);
  //double qdev=qdev_COEF*pow(L12, -qdev_EXP);
  //qinf=qave+A*qdev;
  return(qave);
}

float Compute_Dcont(double *qinf, float q, int L1, int L2,
		    int *related, int NORM)
{
  double L12, qave, qdev, qmax=0.25;
  if(Ini_Dcont==0){
    Ini_Dcont=1; L12=L_max;
    qave=qave_COEF*pow(L12, -qave_EXP);
    qdev=qdev_COEF*pow(L12, -qdev_EXP);
    *qinf=qave; //+A*qdev;
    D0=1.-log(qdev)+log(1.-*qinf);//+A
  }

  if(NORM==0){if(L1<L2){L12=L1;}else{L12=L2;}}
  else if(NORM==1){if(L1>L2){L12=L1;}else{L12=L2;}}
  else{L12=sqrt((float)(L1*L2));}

  qave=qave_COEF*pow(L12, -qave_EXP);
  qdev=qdev_COEF*pow(L12, -qdev_EXP);
  *qinf=qave; //+A*qdev;
  if(*qinf>qmax)*qinf=qmax;
  //printf("qave= %.4f qdev= %.4f qinf= %.4f D0=%.2f L=%.0f\n",
  // qave, qdev, qinf, D0, L12);
  if(q>*qinf+qdev){
    float D1=-log((q-*qinf)/(1.-*qinf));
    *related=1;
    //printf("D1=%.2f D2=%.2f\n", D1, D2);
    //if((z>1)||(D1<D2))
    return(D1);
  }
  float z=(q-qave)/qdev;
  float D2=D0-z;
  *related=0;
  return(D2);
}
