#include "ODESolvers.h"
#include "gnuplot_i.h"
#include "utilities.h"
#include <stdio.h>
#include <string.h>

// for outputfilename
void specifySolverMethodInit(odeOptions *options){

    switch(options -> methodId) {
        case 1: options -> method = "FWE"; break;
        case 2: options -> method = "Heun"; break;
        case 3: options -> method = "Midpoint"; break;
        case 4: options -> method = "RK2Ralston"; break;
        case 5: options -> method = "RK3Classic"; break;
        case 6: options -> method = "RK3Optim"; break;
        case 7: options -> method = "RK4Classic"; break;
        case 8: options -> method = "RK5Butcher"; break;
        default: printf("Incorrect methodId declared. Exiting program..\n"); exit(EXIT_FAILURE);
    }
}

// -- Clear Memory -----------------------------------------------------------

void delete(solution *result, odeOptions *options){

    printf("\n----- MEMORY DEALLOCATION START ->");

    gsl_vector_free(result -> dom);
    gsl_matrix_free(result -> func);
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

    fprintf(outputfile, "#Domain,Functions\n");

    for(largeInt point = 0; point <= options -> lastIndex; ++point){
        fprintf(outputfile, "%012.9lf", gsl_vector_get(result -> dom, point));

        for (int var = 0; var < options -> NSYS; ++var) {
            fprintf(outputfile, ",%012.9lf", gsl_matrix_get(result -> func, var, point));
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
                printf("\t%12.9lf", gsl_matrix_get(result -> func, var, point));
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
