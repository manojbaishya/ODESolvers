#ifndef ODE_SOLVERS_H
#define ODE_SOLVERS_H

#include <stdbool.h>

// -- typedefs and data structures --------------------------------------------

typedef unsigned long long int largeInt;
typedef struct _odeOptions odeOptions;
typedef struct _solution solution;


void callODESolver(void (*)(const double *, const double [], double []), int (*)(const double *, const double []), const char *, int);
odeOptions * readInput(const char *, int);
void ODEinit(odeOptions *, int (*)(const double *, const double []));
void specifySolverMethodInit(odeOptions *);


solution * ODESolver(void (*)(const double *, const double [], double []), odeOptions *);
void ODEIntegrate(void (*)(const double *, const double [], double []), solution *, odeOptions *, largeInt, double);
int adaptiveODEIntegrate(void (*)(const double *, const double [], double []), solution *, odeOptions *, largeInt);

void genericSolver(void (*)(const double *, const double [], double []), double *, double *, double, odeOptions *);

// -- Utility Functions -------------------------------------------------------

void realloc_gsl_containers(solution *, odeOptions *);
void printResult(solution *, odeOptions *);
void writefile(solution *, odeOptions *);
void plotData(solution *, odeOptions *);
void delete(solution *, odeOptions *);

#endif // ODE_SOLVERS_H
