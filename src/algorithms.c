/*
* @Author: manoj
* @Date:   2020-04-03 17:24:45
* @Last Modified by:   Manoj Baishya
* @Last Modified time: 2020-04-26 16:58:03
*/

#include "algorithms.h"
#include <math.h>

// -- Algorithms for first order -------------------------------------------

void FWEuler(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step){
    double slope; derivative(t, y, &slope);
    *y = *y + slope * step;
    *t = *t + step;
}

void Heun(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step){

    double yi = *y; double y_old = 0.0;
    // double phi_avg = 0.0; // = (phi_i + phi_i+1)/2

    double phi_i, phi_ip1;
    derivative(t, y, &phi_i); // slope at i-th point

    *y = yi + phi_i * step; // predictor, y_i+1_0
    *t = *t + step; // t_i+1

    double error = 0.0; double errorBound = 0.5; // in percent, 0.5%
    int iter = 1; int max_iter = 100;

    do {

        y_old = *y;
        derivative(t, y, &phi_ip1);
        *y = yi + ((phi_i + phi_ip1) * step)/2;
        error = fabs((*y - y_old)/(*y)) * 100;
        ++iter;

    } while (error >= errorBound || iter <= max_iter);

}

void Midpoint(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step){

    double slope_st, slope_mid;

    derivative(t, y, &slope_st);

    double t_half = *t + step/2;
    double y_half = *y + slope_st * (step/2);
    derivative(&t_half, &y_half, &slope_mid);

    *y = *y + slope_mid * step;
    *t = *t + step;
}


void RK2Ralston(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step){
    double k1, k2, t_tmp, y_tmp;

    derivative(t, y, &k1);

    t_tmp = *t + 0.75 * step; y_tmp = *y + 0.75 * k1 * step;
    derivative(&t_tmp, &y_tmp, &k2);

    *y = *y + ((1.0/3.0) * k1 + (2.0/3.0) * k2) * step;
    *t = *t + step;
}

void RK3Classic(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step){
    double k1, k2, k3, t_tmp, y_tmp;
    derivative(t, y, &k1);
    t_tmp = *t + 0.5 * step; y_tmp = *y + 0.5 * k1 * step;
    derivative(&t_tmp, &y_tmp, &k2);
    t_tmp = *t + step; y_tmp = *y + (-k1 + 2.0 * k2) * step;
    derivative(&t_tmp, &y_tmp, &k3);

    *y = *y + ((1.0/6.0) * k1 + (2.0/3.0) * k2 + (1.0/6.0) * k3) * step;
    *t = *t + step;
}

void RK3Optim(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step){
    double k1, k2, k3, t_tmp, y_tmp;

    derivative(t, y, &k1);
    t_tmp = *t + (1.0/3.0) * step; y_tmp = *y + (1.0/3.0) * k1 * step;
    derivative(&t_tmp, &y_tmp, &k2);
    t_tmp = *t + (2.0/3.0) * step; y_tmp = *y + (2.0/3.0) * k2 * step;
    derivative(&t_tmp, &y_tmp, &k3);

    *y = *y + ((1.0/4.0) * k1 + (3.0/4.0) * k3) * step;
    *t = *t + step;
}

void RK4Classic(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step){
    double k1, k2, k3, k4, t_tmp, y_tmp;
    derivative(t, y, &k1);
    t_tmp = *t + 0.5 * step; y_tmp = *y + 0.5 * k1 * step;
    derivative(&t_tmp, &y_tmp, &k2);
    y_tmp = *y + 0.5 * k2 * step;
    derivative(&t_tmp, &y_tmp, &k3);
    t_tmp = *t + step; y_tmp = *y + k3 * step;
    derivative(&t_tmp, &y_tmp, &k4);

    *y = *y + (step/6.0) * (k1 + 2.0 * (k2 + k3) + k4);
    *t = *t + step;
}

void RK5Butcher(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double *y, double step){
    double k1, k2, k3, k4, k5, k6, t_tmp, y_tmp;

    derivative(t, y, &k1);
    t_tmp = *t + 0.25 * step; y_tmp = *y + 0.25 * k1 * step;
    derivative(&t_tmp, &y_tmp, &k2);
    t_tmp = *t + 0.25 * step; y_tmp = *y + (0.125 * k1 + 0.125 * k2) * step;
    derivative(&t_tmp, &y_tmp, &k3);
    t_tmp = *t + 0.5 * step; y_tmp = *y + (-0.5 * k2 + k3) * step;
    derivative(&t_tmp, &y_tmp, &k4);
    t_tmp = *t + 0.75 * step; y_tmp = *y + ((3.0/16.0) * k1 + (9.0/16.0) * k4) * step;
    derivative(&t_tmp, &y_tmp, &k5);
    t_tmp = *t + step; y_tmp = *y + (-(3.0/7.0) * k1 + (2.0/7.0) * k2 + (12.0/7.0) * k3 - (12.0/7.0) * k4 + (8.0/7.0) * k5) * step;
    derivative(&t_tmp, &y_tmp, &k6);

    *y = *y + (step/90.0) * (7.0 * k1 + 32.0 * k3 + 12.0 * k4 + 32.0 * k5 + 7.0 * k6);
    *t = *t + step;
}

// ----------------------------------------------------------------------------
//
//
//
// ----------------------------------------------------------------------------

void RK4SYS(void (*derivative)(const double *t, const double y[], double ydot[]), double *t, double y[], double step, int NSYS){

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

// ----------------------------------------------------------------------------
//
//
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
