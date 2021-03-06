/*
* @Author: manoj
* @Date:   2020-04-01 14:04:10
* @Last Modified by:   Manoj Baishya
* @Last Modified time: 2020-05-06 16:32:58
*/

#include "ODESolvers.h"
#include "algorithms.h"
#include "utilities.h"
#include "parson.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

// -- Caller Function ---------------------------------------------------------

void callODESolver(void (*derivative)(const double *t, const double y[], double ydot[]), int (*events)(const double *, const double []), const char *inputfile, int NSYS){

    puts("\n---------------------- Starting the program! ----------------------\n");

    printf("\t- Solving system of ODEs...\n");

    odeOptions *options = readInput(inputfile, NSYS);

    ODEinit(options, events);

    solution *result = ODESolver(derivative, options);

    // post-process data

    writefile(result, options);

    printResult(result, options);

    plotData(result, options);

    delete(result, options);

    // clear memory

    printf("\n---------------------- EXITING PROGRAM ----------------------\n");
}


// -- Input Reader Function ---------------------------------------------------

odeOptions * readInput(const char *inputjson, int NSYS){

    odeOptions *options = (odeOptions *) malloc(sizeof(odeOptions) + sizeof(long double) * NSYS);
    JSON_Value *file = json_parse_file_with_comments(inputjson);
    JSON_Object *data = json_object(file);

    JSON_Array *buffer = json_object_get_array(data, "domain");
    size_t count = json_array_get_count(buffer);
    for (size_t index = 0; index < count; ++index) {
        options -> domain[index] = json_array_get_number(buffer, index);
    }

    options -> step = json_object_get_number(data, "stepsize");
    options -> outInterval = json_object_get_number(data, "outputInterval");
    options -> relErr = json_object_get_number(data, "relative_errorPC") / 100;
    options -> NSYS = NSYS;
    options -> adaptive = (bool) json_object_get_number(data, "adaptive_switch");
    options -> methodId = json_object_get_number(data, "methodId");

    options -> printResult = json_object_get_number(data, "printResult");
    options -> plotTimeSeries = json_object_get_number(data, "plotTimeSeries");

    options -> model = (char *) malloc(sizeof(char) * (strlen(json_object_get_string(data, "modelname")) + 1));
    strcpy(options -> model, json_object_get_string(data, "modelname"));

    buffer = json_object_get_array(data, "yInitCond");
    count = json_array_get_count(buffer);
    for (size_t index = 0; index < count; ++index) {
        options -> yInitCond[index] = json_array_get_number(buffer, index);
    }

    json_value_free(file);

    return options;
}


// -- Initialisation Function -------------------------------------------------

void ODEinit(odeOptions *options, int (*events)(const double *, const double [])){

    // select solver method
    (options -> adaptive == 1) ? (options -> method = "CashKarpRKF45") : specifySolverMethodInit(options);

    options -> GRIDPOINTS = (largeInt) ((options -> domain[1] - options -> domain[0])/options -> outInterval) + 1;

    // specify outputfilepath ----------------

    char filepath[100] = "./workspace/data/";
    strcat(filepath, options -> model);

    // Create subdirectory in ./data if None
    struct stat st = {0};
    if (stat(filepath, &st) == -1) {
       mkdir(filepath, S_IRWXU);
    }

    sprintf(&filepath[strlen(filepath)], "/%s_step=%4.2le.csv", options -> method, options -> step);

    options -> outputFilePath = (char *) malloc(sizeof(char) * (strlen(filepath) + 1));
    strcpy(options -> outputFilePath, filepath);

    // Assign events function pointer
    options -> events = events;

}

// -- Templated Solvers --------------------------------------------------------

solution * ODESolver(void (*derivative)(const double *t, const double y[], double ydot[]), odeOptions *options){

    printf("\n\t- Using %s algorithm!\n\t- Solution in progress...\n\n", options -> method);

    // allocate memory for storing results
    solution *result = (solution *) malloc(sizeof(solution));

    result -> dom = gsl_vector_alloc(options -> GRIDPOINTS);
    result -> func = gsl_matrix_alloc(options -> NSYS, options -> GRIDPOINTS);

    // assign initial conditions
    gsl_vector_set(result -> dom, 0, options -> domain[0]);
    for (int var = 0; var < options -> NSYS; ++var) {
        gsl_matrix_set(result -> func, var, 0, options -> yInitCond[var]);
    }

    double endtime = 0.0;
    largeInt point;
    int eventflag = 0;

    for (point = 0; gsl_vector_get(result -> dom, point) < options -> domain[1]; ++point) {

        if(options -> adaptive == 1) {

            // Realloc memory if array bounds exceeded

            if(point + 1 == options -> GRIDPOINTS) {
                realloc_gsl_containers(result, options);
            }

            eventflag = adaptiveODEIntegrate(derivative, result, options, point);

            if(eventflag == 1) {
                break;
            }

        } else if(options -> adaptive == 0) {

            // Realloc memory if array bounds exceeded

            if(point + 1 == options -> GRIDPOINTS) {
                realloc_gsl_containers(result, options);
            }

            endtime = gsl_vector_get(result -> dom, point) + options -> outInterval;

            if(endtime > options -> domain[1]) {
                endtime = options -> domain[1];
            }

            ODEIntegrate(derivative, result, options, point, endtime);
        }
    }

    options -> lastIndex = point - 1;

    puts("---------------------- ODE solved successfully! ----------------------\n");

    return result;

}

// -- Single Step Integrator --------------------------------------------------


int adaptiveODEIntegrate(void (*derivative)(const double *t, const double y[], double ydot[]), solution *result, odeOptions *options, largeInt point) {

    // extract solution to temporary containers
    double step = options -> step;
    double indep_t = gsl_vector_get(result -> dom, point);
    double func_y[options -> NSYS];
    for (int var = 0; var < options -> NSYS; ++var) {
        func_y[var] = gsl_matrix_get(result -> func, var, point);
    }

    // specify algorithm parameters
    // {'nonZeroScaffold': this parameter guards against driving step size to zero (infinitesimal)! }
    static double safety = 0.9, errorMinBound = 5.7665e-4, nonZeroScaffold = 1.0e-30;
    double errorMax;
    double ytemp[options -> NSYS], errorSpectrum[options -> NSYS];
    double dydt[options -> NSYS], yscal[options -> NSYS];

    derivative(&indep_t, func_y, dydt);
    for (int var = 0; var < options -> NSYS; ++var) {
        // around y[i] == 0, h * dydt = finite and nonZeroScaffold > 0 implies yscal[i] doesn't go to zero!
        // for high values of y[i], h * dydt and nonZeroScaffold are negligible
        yscal[var] = fabs(func_y[var]) + fabs(step * dydt[var]) + nonZeroScaffold;
    }

    if(indep_t + step > options -> domain[1]) {
        step = options -> domain[1] - indep_t;
    }

    while(true) {

        // trial solution
        CashKarp_RKF45(derivative, &indep_t, func_y, ytemp, step, errorSpectrum, options -> NSYS);

        // determine error signal from trial solution
        errorMax = 0.0;
        for (int var = 0; var < options -> NSYS; ++var){
            errorMax = FMAX(errorMax, fabs(errorSpectrum[var]/yscal[var]));
        }
        errorMax /= options -> relErr;

        // step modification based on error feedback
        if(errorMax > 1.0) {

            // scale down stepsize by a maximum factor of 4
            step = FMAX(fabs(safety * step * pow(errorMax, -0.25)), 0.25 * fabs(step));

            if(indep_t + step == indep_t) {
                fprintf(stderr, "\nstepsize underflow in adaptive stepper algorithm..now exiting to system\n");
                exit(1);
            }

            // continue inner loop (repeat trial solution)
        } else {

            // errorMax driven to less than 1

            indep_t = indep_t + step; // advance time step
            for (int var = 0; var < options -> NSYS; ++var) {
                func_y[var] = ytemp[var]; // advance solution by finer time step
            }

            // scale up stepsize by a maximum factor of 4 only
            step = (errorMax > errorMinBound) ? safety * step * pow(errorMax, -0.20) : (4.0 * step); // valid stepsize for next step

            break; // out of inner loop, solution for this step successful based on specified error
        }
    } // end inner while loop

    // ===================== Check Event =========================
    int eventflag = 0;
    double eventTime, eventSol[options -> NSYS];

    eventTime = gsl_vector_get(result -> dom, point);

    for (int var = 0; var < options -> NSYS; ++var) {
        eventSol[var] = (double) gsl_matrix_get(result -> func, var, point);
    }

    eventflag = options -> events(&eventTime, eventSol);

    if(eventflag == 0) {
        options -> step = step; // valid stepsize for next outer coarse loop step
        gsl_vector_set(result -> dom, point + 1, indep_t);

        for (int var = 0; var < options -> NSYS; ++var) {
            gsl_matrix_set(result -> func, var, point + 1, func_y[var]);
        }
    }

    return eventflag;

}


void ODEIntegrate(void (*derivative)(const double *t, const double y[], double ydot[]), solution *result, odeOptions *options, largeInt point, double endtime){

    double step = options -> step;
    double indep_t = gsl_vector_get(result -> dom, point);

    double func_y[options -> NSYS];
    for (int var = 0; var < options -> NSYS; ++var) {
        func_y[var] = gsl_matrix_get(result -> func, var, point);
    }

    do {

        if(endtime - indep_t < step) {
            step = endtime - indep_t;
        }

        genericSolver(derivative, &indep_t, func_y, step, options);

    } while (indep_t < endtime);

    gsl_vector_set(result -> dom, point + 1, indep_t);

    for (int var = 0; var < options -> NSYS; ++var) {
        gsl_matrix_set(result -> func, var, point + 1, func_y[var]);

    }
}

void genericSolver(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step, odeOptions *options){

    switch(options -> methodId) {
        case 1: FWEuler(derivative, t, y, step, options -> NSYS); break;
        case 2: Heun(derivative, t, y, step, options -> NSYS); break;
        case 3: Midpoint(derivative, t, y, step, options -> NSYS); break;
        case 4: RK2Ralston(derivative, t, y, step, options -> NSYS); break;
        case 5: RK3Classic(derivative, t, y, step, options -> NSYS); break;
        case 6: RK3Optim(derivative, t, y, step, options -> NSYS); break;
        case 7: RK4(derivative, t, y, step, options -> NSYS); break;
        case 8: RK5Butcher(derivative, t, y, step, options -> NSYS); break;
    }

}

void realloc_gsl_containers(solution *result, odeOptions *options){

    largeInt GRIDPOINTS = options -> GRIDPOINTS + 50;

    gsl_vector *temp_domain = gsl_vector_alloc(GRIDPOINTS);
    gsl_matrix *temp_func = gsl_matrix_alloc(options -> NSYS, GRIDPOINTS);

    for(largeInt point = 0; point < options -> GRIDPOINTS; ++point){
        gsl_vector_set(temp_domain, point, gsl_vector_get(result -> dom, point));
        for (int var = 0; var < options -> NSYS; ++var) {
            gsl_matrix_set(temp_func, var, point, gsl_matrix_get(result -> func, var, point));
        }
    }

    gsl_vector_free(result -> dom); gsl_matrix_free(result -> func);

    result -> dom = temp_domain;
    result -> func = temp_func;

    options -> GRIDPOINTS = GRIDPOINTS;
}
