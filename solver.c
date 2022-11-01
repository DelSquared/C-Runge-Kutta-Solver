#include <stdio.h>
#include <stdlib.h>

void ode(float* y, float t, float* dy, int N){
  for (int n = 0; n < N; n++){
    dy[n] = -2*y[n];
  }
}

void RungeKuttaStep(float* y0, float t, float* y1, int N, float dt){
  float* _y = malloc(N);
  float* k1 = malloc(N);
  ode(y0, t, k1, N);
  float* k2 = malloc(N);
  for (int n = 0; n < N; n++) _y[n] = y0[n]+k1[n]*dt/2;
  ode(_y, t+dt/2, k2, N);
  float* k3 = malloc(N);
  for (int n = 0; n < N; n++) _y[n] = y0[n]+k2[n]*dt/2;
  ode(_y, t+dt/2, k3, N);
  float* k4 = malloc(N);
  for (int n = 0; n < N; n++) _y[n] = y0[n]+k3[n]*dt;
  ode(_y, t+dt, k4, N);

  for (int n = 0; n < N; n++) y1[n] = y0[n] + (k1[n] + 2*k2[n] + 2*k3[n] + k4[n])*dt/6;

  free(_y);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
}

int main() {
  int N = 4;
  float t = 2;
  float dt = 0.01;
  float* y0 = malloc(N);
  float* y1 = malloc(N);
  for (int n = 0; n < N; n++){
    y0[n] = n+1;
  }
  for (int n = 0; n < N; n++){
    printf("%f\n", y0[n]);
  }
  printf("\n");
  RungeKuttaStep(y0, t, y1, N, dt);
  for (int n = 0; n < N; n++){
    printf("%f\n", y1[n]);
  }
  printf("exiting");
  return 0;
}
