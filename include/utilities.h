#ifndef UTILITIES_H
#define UTILITIES_H

#include "ODESolvers.h"

#define FMAX(x, y) ( x > y ? x : y )
#define deg2rad(ang_deg) ((ang_deg * M_PI) / 180.0)
#define rad2deg(ang_rad) ((ang_rad * 180.0) / M_PI)

void specifySolverMethodInit(odeOptions *);
void printResult(solution *, odeOptions *);
void writefile(solution *, odeOptions *);
void plotData(solution *, odeOptions *);
void delete(solution *, odeOptions *);

#endif // UTILITIES_H
