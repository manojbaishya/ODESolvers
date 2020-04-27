/*
* @Author: manoj
* @Date:   2020-03-22 19:03:32
* @Last Modified by:   Manoj Baishya
* @Last Modified time: 2020-04-27 19:41:56
*/

#include "derivatives.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// -- Inline Functions --------------------------------------------------------

#define deg2rad(ang_deg) ((ang_deg * M_PI) / 180.0)
#define rad2deg(ang_rad) ((ang_rad * 180.0) / M_PI)

// ----------------------------------------------------------------------------
//
//                 Required definitions in derivative.c:
//
//  1) g_NSYS
//  2) params
//  3) set_parameters()
//  4) derivative()
//  5) derivative_internal()
//  6) events()
//
// ----------------------------------------------------------------------------



/*
int g_NSYS = 1;

struct params {
    double g, m, c, a, b, vmax;
};

void derivative(const double *t, const double y[], double ydot[]){
    static struct params consts;

    consts.g = 9.81;
    consts.m = 68.1;
    consts.c = 12.5;
    consts.a = 8.3;
    consts.b = 2.2;
    consts.vmax = 46;

    derivative_internal(t, y, ydot, consts);
}

void derivative_internal(const double *t, const double y[], double ydot[], struct params consts){
    ydot[0] = consts.g - ((consts.c)/(consts.m)) * (y[0] + consts.a * pow(y[0]/consts.vmax, consts.b));
}
*/


// ----------------------------------------------------------------------------
//
//
//
// ----------------------------------------------------------------------------


/*
int g_NSYS = 6; // Order of the system of equations

struct params {
    double g, m1, m2, m3, k1, k2, k3;
} g_consts;

void set_parameters(struct params *constptr){
    int choice;
    printf("Enter a value for k2: [1] 100, [2] 150, [3] 200: \n");
    scanf("%d", &choice); getchar();

    switch(choice) {
        case 1: constptr -> k2 = 100; break;
        case 2: constptr -> k2 = 150; break;
        case 3: constptr -> k2 = 200; break;
        default: printf("Incorrect input entered. Exiting program..\n"); exit(EXIT_FAILURE); break;
    }

    constptr -> g = 9.81; // m/s^2
    constptr -> m1 = 60; // kg
    constptr -> m2 = 70; // kg
    constptr -> m3 = 80; // kg
    constptr -> k1 = g_consts.k3 = 50; // N/m
}

void derivative(const double *t, const double y[], double ydot[]){

    derivative_internal(t, y, ydot, g_consts);

}

void derivative_internal(const double *t, const double y[], double ydot[], struct params consts){

    ydot[0] = y[1];
    ydot[1] = consts.g + (consts.k2 * (y[2] - y[0]) - consts.k1 * y[0])/consts.m1;
    ydot[2] = y[3];
    ydot[3] = consts.g + (consts.k3 * (y[4] - y[2]) + consts.k2 * (y[0] - y[2]))/consts.m2;
    ydot[4] = y[5];
    ydot[5] = consts.g + (consts.k3 * (y[2] - y[4]))/consts.m3;
}
*/

// ----------------------------------------------------------------------------
//
//
//
// ----------------------------------------------------------------------------

/*
int g_NSYS = 6; // Order of the system of equations !Required

struct params {
    double Vt, Vm, AlphaT, del;
} g_consts;

void set_parameters(struct params *constptr){

    int choice = 0;
    printf("\nEnter a value for K: [1] 0.83, [2] 1.11, [3] 1.67, [4] 5: ");
    scanf("%d", &choice); getchar();

    float K = 0.0;
    switch(choice) {
        case 1: K = 0.83; break;
        case 2: K = 1.11; break;
        case 3: K = 1.67; break;
        case 4: K = 5; break;
        default: printf("Incorrect input entered. Exiting program..\n"); exit(EXIT_FAILURE); break;
    }

    constptr -> AlphaT = M_PI;
    constptr -> del = deg2rad(30);
    constptr -> Vt = 300;
    constptr -> Vm = K * constptr -> Vt;
}

void derivative(const double *t, const double y[], double ydot[]){

    derivative_internal(t, y, ydot, g_consts);

}

void derivative_internal(const double *t, const double y[], double ydot[], struct params consts){

    // y[0] = R, y[1] = Theta, y[2] = Xm, y[3] = Ym, y[4] = Xt, y[5] = Yt

    ydot[0] = (consts.Vt) * cos(consts.AlphaT - y[1]) - consts.Vm * cos(consts.del);
    ydot[1] = (consts.Vt * sin(consts.AlphaT - y[1]) - consts.Vm * sin(consts.del))/y[0];
    ydot[2] = consts.Vm * cos(y[1] + consts.del);
    ydot[3] = consts.Vm * sin(y[1] + consts.del);
    ydot[4] = consts.Vt * cos(consts.AlphaT);
    ydot[5] = consts.Vt * sin(consts.AlphaT);
}

int events(const double *t, const double y[]) {

    int numEvents = 1;
    double event_checks[] = {y[0] - 0.001}; // count of values in initialiser must equal numEvents

    int flag = 0;

    for(unsigned i = 0; i < numEvents; ++i) {
        if( event_checks[i] < 0 ) {
            flag = 1;
            break;
        } else {
            flag = 0;
        }
    }

    return flag;
}
*/
// ----------------------------------------------------------------------------
//
//
//
// ----------------------------------------------------------------------------



int g_NSYS = 1;

struct params {
    double k, mu, sig;
} g_consts;


void set_parameters(struct params *constptr){

    constptr -> mu = 2.0;
    constptr -> sig = 0.075;
    constptr -> k = 0.6;
}


void derivative(const double *t, const double y[], double ydot[]){

    derivative_internal(t, y, ydot, g_consts);
}

void derivative_internal(const double *t, const double y[], double ydot[], struct params consts){
    ydot[0] = - consts.k * y[0] + 10 * exp(- (pow((*t - consts.mu), 2)/(2 * pow(consts.sig, 2))));
}

int events(const double *t, const double y[]) {
/**
 * @brief Returns event information
 * @details set all flag returns to 0 if there is no event
 *
 * @param t time
 * @param y[] solution at time t
 * @return eventflag
 */
    return 0;
}



