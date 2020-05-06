#ifndef ALGORITHMS_H
#define ALGORITHMS_H

// -- nth Order Systems, Non-Adaptive ----------------------------------------------------

void FWEuler(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS);
void Heun(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS);
void Midpoint(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS);
void RK2Ralston(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS);
void RK3Classic(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS);
void RK3Optim(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS);
void RK4(void (*)(const double *, const double [], double []), double *, double [], double, int);
void RK5Butcher(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS);

// -- nth Order Systems, Adaptive ----------------------------------------------------

void CashKarp_RKF45(void (*)(const double *, const double [], double []), double *, double [], double [], double, double [], int);

// -- Root Finder ----------------------------------------------------------//
double newton_raphson(double (*)(double), double (*)(double), double, int);

#endif // ALGORITHMS_H
