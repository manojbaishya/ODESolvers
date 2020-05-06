/*
* @Author: manoj
* @Date:   2020-04-03 17:24:45
* @Last Modified by:   Manoj Baishya
* @Last Modified time: 2020-05-06 16:32:22
*/

#include "algorithms.h"
#include <stdio.h>
#include <math.h>

// -- Macro/Inline Functions ---------------------------------------------------------

#define FMAX(x, y) ( x > y ? x : y )

// ----------------------------------------------------------------------------
//
//                            Non Adaptive systems algorithms
//
// ----------------------------------------------------------------------------

// Method ID = 1
void FWEuler(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS){
    double slope[NSYS]; derivative(t, y, slope);
    for (int index = 0; index < NSYS; ++index) {
        y[index] = y[index] + slope[index] * step;
    }
    *t = *t + step;
}

// Method ID = 2
void Heun(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS){

    double yi[NSYS];
    double y_old[NSYS];

    for (int index = 0; index < NSYS; ++index) {
        yi[index] = y[index];
    }


    double phi_i[NSYS], phi_ip1[NSYS];
    derivative(t, y, phi_i); // slope at i-th point

    for (int index = 0; index < NSYS; ++index) {
        y[index] = yi[index] + phi_i[index] * step; // predictor, y_i+1_0
    }
    *t = *t + step; // t_i+1

    double error[NSYS];
    double errorNorm, errorBound = 0.01; // in percent, 0.5%
    int iter = 1, max_iter = 200;

    do {
        for (int index = 0; index < NSYS; ++index) {
            y_old[index] = y[index];
        }

        derivative(t, y, phi_ip1);

        for (int index = 0; index < NSYS; ++index) {
            y[index] = yi[index] + ((phi_i[index] + phi_ip1[index]) * step)/2;
            error[index] = fabs((y[index] - y_old[index])/y[index]) * 100;
        }

        // L-inf Norm
        errorNorm = 0.0;
        for (int ind = 0; ind < NSYS; ++ind) {
            errorNorm = FMAX(errorNorm, error[ind]);
        }

        iter++;

    } while(iter <= max_iter && errorNorm >= errorBound);


}

// Method ID = 3
void Midpoint(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step, int NSYS){

    double slope_st[NSYS], slope_mid[NSYS];
    double t_half, y_half[NSYS];

    derivative(t, y, slope_st);

    t_half = *t + step/2;
    for (int index = 0; index < NSYS; ++index) {
        y_half[index] = y[index] + slope_st[index] * (step/2);
    }

    derivative(&t_half, y_half, slope_mid);

    for (int index = 0; index < NSYS; ++index) {
        y[index] = y[index] + slope_mid[index] * step;
    }

    *t = *t + step;

}

// Method ID = 4
void RK2Ralston(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step, int NSYS){

    double k1[NSYS], k2[NSYS], t_tmp, y_tmp[NSYS];

    derivative(t, y, k1);

    t_tmp = *t + 0.75 * step;

    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + 0.75 * k1[index] * step;
    }

    derivative(&t_tmp, y_tmp, k2);

    for (int index = 0; index < NSYS; ++index) {
        y[index] = y[index] + ((1.0/3.0) * k1[index] + (2.0/3.0) * k2[index]) * step;
    }

    *t = *t + step;
}

// Method ID = 5
void RK3Classic(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step, int NSYS){

    double k1[NSYS], k2[NSYS], k3[NSYS], t_tmp, y_tmp[NSYS];
    derivative(t, y, k1);

    t_tmp = *t + 0.5 * step;

    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + 0.5 * k1[index] * step;
    }
    derivative(&t_tmp, y_tmp, k2);

    t_tmp = *t + step;
    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + (-k1[index] + 2.0 * k2[index]) * step;
    }

    derivative(&t_tmp, y_tmp, k3);

    for (int index = 0; index < NSYS; ++index) {
        y[index] = y[index] + ((1.0/6.0) * k1[index] + (2.0/3.0) * k2[index] + (1.0/6.0) * k3[index]) * step;
    }

    *t = *t + step;
}


// Method ID = 6
void RK3Optim(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step, int NSYS){

    double k1[NSYS], k2[NSYS], k3[NSYS], t_tmp, y_tmp[NSYS];

    derivative(t, y, k1);

    t_tmp = *t + (1.0/3.0) * step;
    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + (1.0/3.0) * k1[index] * step;
    }
    derivative(&t_tmp, y_tmp, k2);

    t_tmp = *t + (2.0/3.0) * step;
    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + (2.0/3.0) * k2[index] * step;
    }

    derivative(&t_tmp, y_tmp, k3);

    for (int index = 0; index < NSYS; ++index) {
        y[index] = y[index] + ((1.0/4.0) * k1[index] + (3.0/4.0) * k3[index]) * step;
    }
    *t = *t + step;
}


// Method ID = 7
void RK4(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS){

    double t_int, y_int[NSYS];

    // Block 1 Calculations -------------------
    double K1[NSYS]; derivative(t, y, K1);


    // Block 2 Calculations -------------------
    t_int = *t + (step/2.0);
    for (int index = 0; index < NSYS; ++index) {
        y_int[index] = y[index] + (step/2.0) * K1[index];
    }

    double K2[NSYS]; derivative(&t_int, y_int, K2);


    // Block 3 Calculations -------------------
    t_int = *t + (step/2.0);
    for (int index = 0; index < NSYS; ++index) {
        y_int[index] = y[index] + (step/2.0) * K2[index];
    }

    double K3[NSYS]; derivative(&t_int, y_int, K3);


    // Block 4 Calculations -------------------
    t_int = *t + step;
    for (int index = 0; index < NSYS; ++index) {
        y_int[index] = y[index] + step * K3[index];
    }

    double K4[NSYS]; derivative(&t_int, y_int, K4);


    // i+1 Increment Step
    double slope[NSYS];

    for (int index = 0; index < NSYS; ++index) {
        slope[index] = (1.0/6.0) * (K1[index] + 2.0 * (K2[index] + K3[index]) + K4[index]);
        y[index] = y[index] + slope[index] * step;
    }

    *t = *t + step;
}


// Method ID = 8
void RK5Butcher(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step, int NSYS){

    double k1[NSYS], k2[NSYS], k3[NSYS], k4[NSYS], k5[NSYS], k6[NSYS], t_tmp, y_tmp[NSYS];

    derivative(t, y, k1);

    t_tmp = *t + 0.25 * step;
    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + 0.25 * k1[index] * step;
    }

    derivative(&t_tmp, y_tmp, k2);

    t_tmp = *t + 0.25 * step;

    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + (0.125 * k1[index] + 0.125 * k2[index]) * step;
    }

    derivative(&t_tmp, y_tmp, k3);

    t_tmp = *t + 0.5 * step;

    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + (-0.5 * k2[index] + k3[index]) * step;
    }

    derivative(&t_tmp, y_tmp, k4);

    t_tmp = *t + 0.75 * step;

    for (int index = 0; index < NSYS; ++index) {
        y_tmp[index] = y[index] + ((3.0/16.0) * k1[index] + (9.0/16.0) * k4[index]) * step;
    }

    derivative(&t_tmp, y_tmp, k5);

    t_tmp = *t + step;

    for (int index = 0; index < NSYS; ++index) {

        y_tmp[index] = y[index] + (-(3.0/7.0) * k1[index] + (2.0/7.0) * k2[index] + (12.0/7.0) * k3[index] - (12.0/7.0) * k4[index] + (8.0/7.0) * k5[index]) * step;
    }

    derivative(&t_tmp, y_tmp, k6);

    for (int index = 0; index < NSYS; ++index) {
        y[index] = y[index] + (step/90.0) * (7.0 * k1[index] + 32.0 * k3[index] + 12.0 * k4[index] + 32.0 * k5[index] + 7.0 * k6[index]);
    }
    *t = *t + step;
}


// ----------------------------------------------------------------------------
//
//                            Adaptive algorithms
//
// ----------------------------------------------------------------------------

void CashKarp_RKF45(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double ytemp[], double step, double errorSpectrum[], int NSYS){

    // -- Parameters ----------------------------------------------------------
    static double
        p1 = 0.2, q11 = 0.2,
        p2 = 0.3, q21 = 3.0/40.0, q22 = 9.0/40.0,
        p3 = 0.6, q31 = 0.3, q32 = -0.9, q33  = 1.2,
        p4 = 1.0, q41 = -11.0/54.0, q42 = 2.5, q43 = -70.0/27.0, q44 = 35.0/27.0,
        p5 = 0.875, q51 = 1631.0/55295.0, q52 = 175.0/512.0, q53 = 575.0/13824.0, q54 = 44275.0/110592.0, q55 = 253.0/4096.0,
        a1 = 37.0/378.0, a3 = 250.0/621.0, a4 = 125.0/594.0, a6 = 512.0/1771.0,
        d1 = (37.0/378.0) - (2825.0/27648.0), d3 = (250.0/621.0) - (18575.0/48384.0),
        d4 = (125.0/594.0) - (13525.0/55296.0), d5 = -(277.0/14336.0), d6 = (512.0/1771.0) - 0.25;

    // --  ------------------------------------------------------------------------

    double t_int, y_int[NSYS];

    // Block 1 Calculations -------------------
    double K1[NSYS]; derivative(t, y, K1);

    // Block 2 Calculations -------------------
    t_int = *t + p1 * step;
    for (int index = 0; index < NSYS; ++index) {
        y_int[index] = y[index] + (q11 * K1[index]) * step;
    }

    double K2[NSYS]; derivative(&t_int, y_int, K2);

    // Block 3 Calculations -------------------
    t_int = *t + p2 * step;
    for (int index = 0; index < NSYS; ++index) {
        y_int[index] = y[index] + (q21 * K1[index] + q22 * K2[index]) * step;
    }

    double K3[NSYS]; derivative(&t_int, y_int, K3);

    // Block 4 Calculations -------------------
    t_int = *t + p3 * step;
    for (int index = 0; index < NSYS; ++index) {
        y_int[index] = y[index] + (q31 * K1[index] + q32 * K2[index] + q33 * K3[index]) * step;
    }

    double K4[NSYS]; derivative(&t_int, y_int, K4);

    // Block 5 Calculations -------------------
    t_int = *t + p4 * step;
    for (int index = 0; index < NSYS; ++index) {
        y_int[index] = y[index] + (q41 * K1[index] + q42 * K2[index] + q43 * K3[index] + q44 * K4[index]) * step;
    }

    double K5[NSYS]; derivative(&t_int, y_int, K5);

    // Block 6 Calculations -------------------
    t_int = *t + p5 * step;
    for (int index = 0; index < NSYS; ++index) {
        y_int[index] = y[index] + (q51 * K1[index] + q52 * K2[index] + q53 * K3[index] + q54 * K4[index] + q55 * K5[index]) * step;
    }

    double K6[NSYS]; derivative(&t_int, y_int, K6);

    // i+1 Increment Step and errors
    double slope[NSYS];

    for (int index = 0; index < NSYS; ++index) {
        slope[index] = a1 * K1[index] + a3 * K3[index] + a4 * K4[index] + a6 * K6[index];

        ytemp[index] = y[index] + slope[index] * step;

        errorSpectrum[index] = (d1 * K1[index] + d3 * K3[index] + d4 * K4[index] + d5 * K5[index] + d6 * K6[index]) * step;
    }

}

// ----------------------------------------------------------------------------
//
//                            Miscellaneous algorithms
//
// ----------------------------------------------------------------------------

// -- Single Root Finder ------------------------------------------------------

// Definitions required:
// double error_func(double abscissa);
// double error_derivative(double abscissa);

double newton_raphson(double (*error_func)(double), double (*error_derivative)(double), double init_guess, int MAX_ITER){

    int count = 1;
    double new_guess, error, root;

    do{
        printf("\nIteration : %d | Root: %3.4f", count, init_guess);

        new_guess = init_guess - error_func(init_guess)/error_derivative(init_guess);
        error = fabs(new_guess - init_guess);

        if(count > MAX_ITER){
            printf("%d iteration counts reached. Exiting loop.\n", MAX_ITER);
            break;
        } else {
            init_guess = new_guess;
            count++;
        }

    } while(error > pow(10, -8));

    printf("\n\nRoot has converged!\n");
    root = new_guess;
    return root;
}
