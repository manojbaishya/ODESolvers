/*
* @Author: manoj
* @Date:   2020-04-01 14:04:10
* @Last Modified by:   Manoj Baishya
* @Last Modified time: 2020-04-27 20:06:14
*/

#include "ODESolvers.h"
#include "algorithms.h"
#include "derivatives.h"
#include "gnuplot_i.h"
#include "parson.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix_long_double.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

// -- Macro/Inline Functions ---------------------------------------------------------

#define FMAX(x, y) ( x > y ? x : y )

// -- Data Structures ---------------------------------------------------------

struct _solution {
    gsl_vector *dom;
    gsl_matrix_long_double *func;
};

struct _odeOptions {
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
    char *method;
    char *outputFilePath;
    double domain[2];
    long double yInitCond[];
};

// ----------------------------------------------------------------------------
//
//                           Solver Setup Functions
//
// ----------------------------------------------------------------------------

// -- Caller Function ---------------------------------------------------------

void callODESolver(void (*derivative)(const double *t, const double y[], double ydot[]), const char *inputfile, int NSYS){

    puts("\n---------------------- Starting the program! ----------------------\n");

    printf("\t- Solving system of ODEs...\n");

    odeOptions *options = readInput(inputfile, NSYS);

    ODEinit(options);

    solution *result = ODESolver(derivative, options);

    writefile(result, options);

    printResult(result, options);

    plotData(result, options);

    delete(result, options);

    printf("\n---------------------- EXITING PROGRAM ----------------------\n");
}


// -- Input Reader Function ---------------------------------------------------

odeOptions * readInput(const char *inputjson, int NSYS){

    odeOptions *options = (odeOptions *) malloc(sizeof(odeOptions) + sizeof(long double) * NSYS);
    JSON_Value *file = json_parse_file_with_comments(inputjson);
    JSON_Object *data = json_object(file);

    // ===========================================================
    JSON_Array *buffer = json_object_get_array(data, "domain");
    size_t count = json_array_get_count(buffer);
    for (size_t index = 0; index < count; ++index) {
        options -> domain[index] = json_array_get_number(buffer, index);
    }
    // =================================================================

    options -> step = json_object_get_number(data, "stepsize");
    options -> outInterval = json_object_get_number(data, "outputInterval");
    options -> relErr = json_object_get_number(data, "relative_errorPC") / 100;
    options -> NSYS = NSYS;
    options -> adaptive = (bool) json_object_get_number(data, "adaptive_switch");

    options -> printResult = json_object_get_number(data, "printResult");
    options -> plotTimeSeries = json_object_get_number(data, "plotTimeSeries");

    options -> model = (char *) malloc(sizeof(char) * (strlen(json_object_get_string(data, "modelname")) + 1));
    strcpy(options -> model, json_object_get_string(data, "modelname"));

    // =================================================================
    buffer = json_object_get_array(data, "yInitCond");
    count = json_array_get_count(buffer);
    for (size_t index = 0; index < count; ++index) {
        options -> yInitCond[index] = json_array_get_number(buffer, index);
    }
    // =================================================================

    json_value_free(file);

    return options;
}


// -- Initialisation Function -------------------------------------------------

void ODEinit(odeOptions *options){

    // select solver method
    (options -> adaptive == 1) ? (options -> method = "CashKarpRKF45") : (options -> NSYS == 1) ? specifySolverMethodInit(options) : (options -> method = "RK4SYS");

    options -> GRIDPOINTS = (largeInt) ((options -> domain[1] - options -> domain[0])/options -> outInterval) + 1;

    // specify outputfilepath ----------------

    char filepath[100] = "./iodata/";
    strcat(filepath, options -> model);

    // Create subdirectory in ./iodata if None
    struct stat st = {0};
    if (stat(filepath, &st) == -1) {
       mkdir(filepath, S_IRWXU);
    }

    sprintf(&filepath[strlen(filepath)], "/%s_step=%4.2le.dat", options -> method, options -> step);

    options -> outputFilePath = (char *) malloc(sizeof(char) * (strlen(filepath) + 1));
    strcpy(options -> outputFilePath, filepath);

}

// ----------------------------------------------------------------------------
//
//                             Interface Functions
//
// ----------------------------------------------------------------------------

// -- Templated Solvers --------------------------------------------------------

solution * ODESolver(void (*derivative)(const double *t, const double y[], double ydot[]), odeOptions *options){

    printf("\n\t- Using %s algorithm!\n\t- Solution in progress...\n\n", options -> method);

    // allocate memory for storing results
    solution *result = (solution *) malloc(sizeof(solution));

    result -> dom = gsl_vector_alloc(options -> GRIDPOINTS);
    result -> func = gsl_matrix_long_double_alloc(options -> NSYS, options -> GRIDPOINTS);

    // assign initial conditions
    gsl_vector_set(result -> dom, 0, options -> domain[0]);
    for (int var = 0; var < options -> NSYS; ++var) {
        gsl_matrix_long_double_set(result -> func, var, 0, options -> yInitCond[var]);
    }

    double tEnd = 0.0;
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
                point = point - 1;
                break;
            }

        } else if(options -> adaptive == 0) {

            tEnd = gsl_vector_get(result -> dom, point) + options -> outInterval;

            if(tEnd > options -> domain[1]) {
                tEnd = options -> domain[1];
            }

            ODEIntegrate(derivative, result, options, point, tEnd);
        }
    }

    options -> lastIndex = point;

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
        func_y[var] = gsl_matrix_long_double_get(result -> func, var, point);
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
        eventSol[var] = (double) gsl_matrix_long_double_get(result -> func, var, point);
    }

    eventflag = events(&eventTime, eventSol);

    if(eventflag == 0) {
        options -> step = step; // valid stepsize for next outer coarse loop step
        gsl_vector_set(result -> dom, point + 1, indep_t);

        for (int var = 0; var < options -> NSYS; ++var) {
            gsl_matrix_long_double_set(result -> func, var, point + 1, func_y[var]);
        }
    }

    return eventflag;

}


void ODEIntegrate(void (*derivative)(const double *t, const double y[], double ydot[]), solution *result, odeOptions *options, largeInt point, double tEnd){

    double step = options -> step;
    double indep_t = gsl_vector_get(result -> dom, point);

    double func_y[options -> NSYS];
    for (int var = 0; var < options -> NSYS; ++var) {
        func_y[var] = gsl_matrix_long_double_get(result -> func, var, point);
    }

    do {

        if(tEnd - indep_t < step) {
            step = tEnd - indep_t;
        }

        (options -> NSYS == 1) ? OneDimAlgorithm(derivative, &indep_t, func_y, step, options): RK4SYS(derivative, &indep_t, func_y, step, options -> NSYS);

    } while (indep_t < tEnd);

    gsl_vector_set(result -> dom, point + 1, indep_t);

    for (int var = 0; var < options -> NSYS; ++var) {
        gsl_matrix_long_double_set(result -> func, var, point + 1, func_y[var]);

    }
}

// -- Solver Selectors ---------------------------------------------------------

void specifySolverMethodInit(odeOptions *options){

    int choice = 0;
    printf("\nSelect ODE Solver Method:\n[1] = Explicit Euler Forward\n[2] = Heun Iterative Predictor Corrector\n[3] = Midpoint\n[4] = RK2Ralston\n[5] = RK3Classic\n[6] = RK3Optim\n[7] = RK4Classic\n[8] = RK5Butcher\nAnswer: ");
    scanf("%d", &choice); getchar();

    switch(choice) {
        case 1: options -> method = "FWE"; break;
        case 2: options -> method = "Heun"; break;
        case 3: options -> method = "Midpoint"; break;
        case 4: options -> method = "RK2Ralston"; break;
        case 5: options -> method = "RK3Classic"; break;
        case 6: options -> method = "RK3Optim"; break;
        case 7: options -> method = "RK4Classic"; break;
        case 8: options -> method = "RK5Butcher"; break;
        default: printf("Incorrect input entered. Exiting program..\n"); exit(EXIT_SUCCESS);
    }
}

void OneDimAlgorithm(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step, odeOptions *options){

    if(strcmp(options -> method, "FWE") == 0) {
        FWEuler(derivative, t, y, step);
    } else if(strcmp(options -> method, "Heun") == 0) {
        Heun(derivative, t, y, step);
    } else if(strcmp(options -> method, "Midpoint") == 0) {
        Midpoint(derivative, t, y, step);
    } else if(strcmp(options -> method, "RK2Ralston") == 0) {
        RK2Ralston(derivative, t, y, step);
    } else if(strcmp(options -> method, "RK3Classic") == 0) {
        RK3Classic(derivative, t, y, step);
    } else if(strcmp(options -> method, "RK3Optim") == 0) {
        RK3Optim(derivative, t, y, step);
    } else if(strcmp(options -> method, "RK4Classic") == 0) {
        RK4Classic(derivative, t, y, step);
    } else if(strcmp(options -> method, "RK5Butcher") == 0) {
        RK5Butcher(derivative, t, y, step);
    }
}


// ----------------------------------------------------------------------------
//
//                              Utility Functions
//
// ----------------------------------------------------------------------------


void realloc_gsl_containers(solution *result, odeOptions *options){

    largeInt GRIDPOINTS = options -> GRIDPOINTS + 50;

    gsl_vector *temp_domain = gsl_vector_alloc(GRIDPOINTS);
    gsl_matrix_long_double *temp_func = gsl_matrix_long_double_alloc(options -> NSYS, GRIDPOINTS);

    for(largeInt point = 0; point < options -> GRIDPOINTS; ++point){
        gsl_vector_set(temp_domain, point, gsl_vector_get(result -> dom, point));
        for (int var = 0; var < options -> NSYS; ++var) {
            gsl_matrix_long_double_set(temp_func, var, point, gsl_matrix_long_double_get(result -> func, var, point));
        }
    }

    gsl_vector_free(result -> dom); gsl_matrix_long_double_free(result -> func);

    result -> dom = temp_domain;
    result -> func = temp_func;

    options -> GRIDPOINTS = GRIDPOINTS;
}


// -- Clear Memory -----------------------------------------------------------

void delete(solution *result, odeOptions *options){

    printf("\n----- MEMORY DEALLOCATION START ->");

    gsl_vector_free(result -> dom);
    gsl_matrix_long_double_free(result -> func);
    free(result);

    free(options -> model);
    free(options -> outputFilePath);
    free(options);

    printf(" MEMORY DEALLOCATION COMPLETE ------\n");

}

// -- Output Functions -----------------------------------------------

void writefile(solution *result, odeOptions *options){

    FILE *outputfile = fopen(options -> outputFilePath, "w+");
    if(outputfile == NULL) {
        perror("Couldn't open file. Exiting program...");
        exit(EXIT_FAILURE);
    }

    fprintf(outputfile, "# - Domain ------- Functions --\n");

    for(largeInt point = 0; point < options -> lastIndex; ++point){
        fprintf(outputfile, "%12.9lf", gsl_vector_get(result -> dom, point));

        for (int var = 0; var < options -> NSYS; ++var) {
            fprintf(outputfile, "\t%12.9Lf", gsl_matrix_long_double_get(result -> func, var, point));
        }

        fprintf(outputfile, "\n");
    }

    fclose(outputfile);

    printf("\t- Data written to %s successfully.\n\t- Please use the gnuplot scripts in ./nbscripts/ to plot.\n", options -> outputFilePath);
}

// -- Data Display Functions --------------------------------------------------

void printResult(solution *result, odeOptions *options){

    if(options -> printResult == 1) {

        printf("\n");

        for (largeInt point = 0; point < options -> lastIndex; ++point) {
            printf("%12.9lf", gsl_vector_get(result -> dom, point));
            for (int var = 0; var < options -> NSYS; ++var) {
                printf("\t%12.9Lf", gsl_matrix_long_double_get(result -> func, var, point));
            }
            printf("\n");
        }

        printf("\n");

    } else {

        printf("\n\t- Skipping printing..\n");
    }
}

void plotData(solution *result, odeOptions *options){

    if(options -> plotTimeSeries == 1) {

        FILE *gnuplotrc = fopen ("./include/gnuplotrc", "rb");

        if (gnuplotrc != NULL) {

            fseek(gnuplotrc, 0, SEEK_END);
            long int length = ftell(gnuplotrc);
            fseek(gnuplotrc, 0, SEEK_SET);

            char *plotSettings = (char *) malloc(length + 1);

            if (plotSettings != NULL) {

                fread(plotSettings, 1, length, gnuplotrc);
                plotSettings[length] = '\0';
                fclose(gnuplotrc);

                printf("\nPlotting data..\t");

                char commands[3000];
                sprintf(commands, "plot ");

                for (int var = 1; var <= options -> NSYS; ++var) {
                    (var < options -> NSYS) ?
                    sprintf(&commands[strlen(commands)],  "\"%s\" using 1:%d with linespoints ls %d, ", options -> outputFilePath, var + 1, var):
                    sprintf(&commands[strlen(commands)],  "\"%s\" using 1:%d with linespoints ls %d", options -> outputFilePath, var + 1, var);
                }

                // Begin plotting
                gnuplot_ctrl *plotAxes = gnuplot_init();

                gnuplot_cmd(plotAxes, "%s", plotSettings);
                gnuplot_cmd(plotAxes, commands);

                gnuplot_close(plotAxes);
                // End plotting

                free(plotSettings);

                printf("Done plotting!\n");
            }
        }
    } else {
        printf("\n\t- Skipping plotting..\n");
    }
}
