#ifndef ODE_SOLVERS_H
#define ODE_SOLVERS_H

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// -- typedefs and data structures --------------------------------------------

typedef unsigned long long int largeInt;

typedef struct _solution {
    gsl_vector *dom;
    gsl_matrix *func;
} solution;

typedef struct _odeOptions {
    double step;
    largeInt GRIDPOINTS;
    largeInt lastIndex;
    double outInterval; // in terms of steps
    double relErr; // error tolerance
    bool adaptive; // adaptive algorithm switch
    int NSYS;
    int printResult;
    int plotTimeSeries;
    char *model;
    int methodId;
    char *method;
    char *outputFilePath;
    int (*events)(const double *t, const double y[]);
    double domain[2];
    double yInitCond[];
} odeOptions;


// -- functions --

odeOptions * readInput(const char *, int);
void callODESolver(void (*)(const double *, const double [], double []), int (*)(const double *, const double []), const char *, int);
void ODEinit(odeOptions *, int (*)(const double *, const double []));


solution * ODESolver(void (*)(const double *, const double [], double []), odeOptions *);
int adaptiveODEIntegrate(void (*)(const double *, const double [], double []), solution *, odeOptions *, largeInt);
void ODEIntegrate(void (*)(const double *, const double [], double []), solution *, odeOptions *, largeInt, double);

void genericSolver(void (*)(const double *, const double [], double []), double *, double *, double, odeOptions *);
void realloc_gsl_containers(solution *, odeOptions *);

#endif // ODE_SOLVERS_H
