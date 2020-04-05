#ifndef ALGORITHMS_H
#define ALGORITHMS_H

// -- First Order Systems, Non-Adaptive ----------------------------------------------------

void FWEuler(void (*)(const double *, const double [], double []), double *, double *, double);
void Heun(void (*)(const double *, const double [], double []), double *, double *, double);
void Midpoint(void (*)(const double *, const double [], double []), double *, double *, double);
void RK2Ralston(void (*)(const double *, const double [], double []), double *, double *, double);
void RK3Classic(void (*)(const double *, const double [], double []), double *, double *, double);
void RK3Optim(void (*)(const double *, const double [], double []), double *, double *, double);
void RK4Classic(void (*)(const double *, const double [], double []), double *, double *, double);
void RK5Butcher(void (*)(const double *, const double [], double []), double *, double *, double);

// -- nth Order Systems, Non-Adaptive ----------------------------------------------------


void RK4SYS(void (*)(const double *, const double [], double []), double *, double [], double, int);

// -- nth Order Systems, Adaptive ----------------------------------------------------

void CashKarp_RKF45(void (*)(const double *, const double [], double []), double *, double [], double [], double, double [], int);

#endif // ALGORITHMS_H
