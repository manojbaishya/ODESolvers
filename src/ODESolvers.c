/*
* @Author: manoj
* @Date:   2020-04-01 14:04:10
* @Last Modified by:   manoj
* @Last Modified time: 2020-04-04 23:46:05
*/

#include "ODESolvers.h"
#include "algorithms.h"
#include "gnuplot_i.h"
#include "parson.h"
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
    double *dom;
    long double **func;
};

struct _odeOptions {
    double step;
    largeInt GRIDPOINTS;
    double outInterval; // in terms of steps
    double relErr; // error tolerance
    bool adaptive; // adaptive algorithm switch
    int NSYS;
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

    printf("\nSolving system of ODEs...\n");

    odeOptions *options = readInput(inputfile, NSYS);

    char flag = 'Y';
    while(flag == 'Y' || flag == 'y') {

        ODEinit(options);

        solution *result = ODESolver(derivative, options);

        printResult(result, options);

        writefile(result, options);

        plotData(result, options);

        delete(result, options);

        printf("\nGenerate another dataset? [Y] or [y] = Yes, [other key] = No\nAnswer: ");
        scanf("%c", &flag);

    }

    printf("\nAll data have been written to ./iodata/ directory. Please use the gnuplot script ./nbscripts/visualise.gpi to plot.\n\n");
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
//                             Interface Functions
//
// ----------------------------------------------------------------------------

// -- Templated Solvers --------------------------------------------------------

solution * ODESolver(void (*derivative)(const double *t, const double y[], double ydot[]), odeOptions *options){

    printf("\nSolving differential equations with %s algorithm.\n\n", options -> method);

    // allocate memory for storing results
    solution *result = (solution *) malloc(sizeof(solution));
    int padding = 5;

    result -> dom = (double *) malloc(sizeof(double) * (options -> GRIDPOINTS + padding));
    result -> func = (long double **) malloc(sizeof(long double *) * options -> NSYS);
    for (int var = 0; var < options -> NSYS; ++var) {
        result -> func[var] = malloc(sizeof(long double) * (options -> GRIDPOINTS + padding));
    }

    // assign initial conditions
    result -> dom[0] = options -> domain[0];
    for (int var = 0; var < options -> NSYS; ++var) {
        result -> func[var][0] = options -> yInitCond[var];
    }

    double tEnd = 0.0;

    for (largeInt point = 0; result -> dom[point] < options -> domain[1]; ++point) {

        tEnd = result -> dom[point] + options -> outInterval;

        if(tEnd > options -> domain[1]) {
            tEnd = options -> domain[1];
        }

        (options -> adaptive == 1) ?
        adaptiveODEIntegrate(derivative, result, options, point, tEnd):
        ODEIntegrate(derivative, result, options, point, tEnd);

    }

    puts("---------------------- ODE solved successfully! ----------------------\n");

    return result;

}

// -- Single Step Integrator --------------------------------------------------

void ODEIntegrate(void (*derivative)(const double *t, const double y[], double ydot[]), solution *result, odeOptions *options, largeInt point, double tEnd){

    double step = options -> step;
    double indep_t = result -> dom[point];
    double func_y[options -> NSYS];
    for (int var = 0; var < options -> NSYS; ++var) {
        func_y[var] = result -> func[var][point];
    }

    do {

        if(tEnd - indep_t < step) {
            step = tEnd - indep_t;
        }

        (options -> NSYS == 1) ? OneDimAlgorithm(derivative, &indep_t, func_y, step, options): RK4SYS(derivative, &indep_t, func_y, step, options -> NSYS);

    } while (indep_t < tEnd);

    result -> dom[point + 1] = indep_t;
    for (int var = 0; var < options -> NSYS; ++var) {
        result -> func[var][point + 1] = func_y[var];
    }
}

void adaptiveODEIntegrate(void (*derivative)(const double *t, const double y[], double ydot[]), solution *result, odeOptions *options, largeInt point, double tEnd) {

    // extract solution to temporary containers
    double step = options -> step;
    double indep_t = result -> dom[point];
    double func_y[options -> NSYS];
    for (int var = 0; var < options -> NSYS; ++var) {
        func_y[var] = result -> func[var][point];
    }

    // specify algorithm parameters
    // {'nonZeroScaffold': this parameter guards against driving step size to zero (infinitesimal)! }
    static double safety = 0.9, errorMinBound = 5.7665e-4, nonZeroScaffold = 1.0e-30, errorMax;
    double ytemp[options -> NSYS], errorSpectrum[options -> NSYS],
    dydt[options -> NSYS], yscal[options -> NSYS];


    do {

        derivative(&indep_t, func_y, dydt);
        for (int var = 0; var < options -> NSYS; ++var) {
            // around y[i] == 0, h * dydt = finite and nonZeroScaffold > 0 implies yscal[i] doesn't go to zero!
            // for high values of y[i], h * dydt and nonZeroScaffold are negligible
            yscal[var] = fabs(func_y[var]) + fabs(step * dydt[var]) + nonZeroScaffold;
        }

        // step overshoot guard at specified output points, only decreases size (safe)!
        if(tEnd - indep_t < step) {
            step = tEnd - indep_t;
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

                indep_t = indep_t + step; // advance finer time step
                for (int var = 0; var < options -> NSYS; ++var) {
                    func_y[var] = ytemp[var]; // advance solution by finer time step
                }

                // scale up stepsize by a maximum factor of 4 only
                step = (errorMax > errorMinBound) ? safety * step * pow(errorMax, -0.20) : (4.0 * step); // valid stepsize for next finer step

                break; // out of inner loop, solution for this step successful based on specified error
            }
        } // end inner while loop for finer step progression

    } while (indep_t < tEnd); // end outer do-while loop for coarser step progression


    options -> step = step; // valid stepsize for next outer coarse loop step
    result -> dom[point + 1] = indep_t; // advance coarser time step
    for (int var = 0; var < options -> NSYS; ++var) {
        result -> func[var][point + 1] = func_y[var]; // advance solution by coarser time step
    }
}

// ----------------------------------------------------------------------------
//
//                              Utility Functions
//
// ----------------------------------------------------------------------------

// -- Clear Memory -----------------------------------------------------------

void delete(solution *result, odeOptions *options){

    printf("\n----- MEMORY DEALLOCATION START ->");

    free(result -> dom);

    for (int var = 0; var < options -> NSYS; ++var) {
        free(result -> func[var]);
    }
    free(result -> func);
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

    fprintf(outputfile, "# - Domain ------- Functions -- \n");

    for(largeInt point = 0; point < options -> GRIDPOINTS; ++point){
        fprintf(outputfile, "%12.9lf", result -> dom[point]);

        for (int var = 0; var < options -> NSYS; ++var) {
            fprintf(outputfile, "\t%12.9Lf", result -> func[var][point]);
        }

        fprintf(outputfile, "\n");
    }

    fclose(outputfile);

    printf("\nData written to %s successfully.\n", options -> outputFilePath);
}

// -- Data Display Functions --------------------------------------------------

void printResult(solution *result, odeOptions *options){

    char print_response = 'y';
    printf("Print results? Answer [(y)es or (n)o]: ");
    scanf("%c", &print_response); getchar();

    if(print_response == 'Y' || print_response == 'y') {

        printf("\n");

        for (largeInt point = 0; point < options -> GRIDPOINTS; ++point) {
            printf("%12.9lf", result -> dom[point]);
            for (int var = 0; var < options -> NSYS; ++var) {
                printf("\t%12.9Lf", result -> func[var][point]);
            }
            printf("\n");
        }

        printf("\n");

    } else {

        printf("\nSkipping printing..\n");
    }
}

void plotData(solution *result, odeOptions *options){
    char plot_response = 'y';
    printf("\nPlot results? Answer [(y)es or (n)o]: ");
    scanf("%c", &plot_response); getchar();

    if(plot_response == 'Y' || plot_response == 'y') {
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
        printf("\nSkipping plotting..\n");
    }
}
