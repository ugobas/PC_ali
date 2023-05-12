#include "coordinates.h"
#include <stdio.h>

float dist2(float *x1, float *x2){
  float dx=x1[0]-x2[0], dy=x1[1]-x2[1], dz=x1[2]-x2[2];
  return(dx*dx+dy*dy+dz*dz);
}

void Average_coord(float *x_ave, float **x, int N){
  x_ave[0]=0; x_ave[1]=0; x_ave[2]=0; 
  for(int i=0; i<N; i++){
    x_ave[0]+=x[i][0]; x_ave[1]+=x[i][1]; x_ave[2]+=x[i][2];
  }
  x_ave[0]/=N; x_ave[1]/=N; x_ave[2]/=N; 
}

void Shift_coord(float **x, float *x_ave, int N){
  for(int i=0; i<N; i++){
    x[i][0]-=x_ave[0]; x[i][1]-=x_ave[1]; x[i][2]-=x_ave[2];
  }
}

void Rotate_coord(float **x, int N, float *r){
  for(int i=0; i<N; i++){
    float *xx=x[i], y[3];
    y[0]=r[0]*xx[0]+r[3]*xx[1]+r[6]*xx[2];
    y[1]=r[1]*xx[0]+r[4]*xx[1]+r[7]*xx[2];
    y[2]=r[2]*xx[0]+r[5]*xx[1]+r[8]*xx[2];
    xx[0]=y[0]; xx[1]=y[1]; xx[2]=y[2];
  }
}

void Copy_coordinates(float **xca, float **xca_store, int nres){
  for(int i=0; i<nres; i++){
    xca[i][0]=xca_store[i][0];
    xca[i][1]=xca_store[i][1];
    xca[i][2]=xca_store[i][2];
  }
}
