/*
* @Author: manoj
* @Date:   2020-03-23 12:49:12
* @Last Modified by:   manoj
* @Last Modified time: 2020-04-05 02:24:56
*/
#include "simulations.h"
#include "ODESolvers.h"
#include "derivatives.h"


void singleODE(void){

    set_parameters(&g_consts);
    callODESolver(derivative, "./include/initial-conditions.json", g_NSYS);

}

void systemODE(void){

    set_parameters(&g_consts);
    callODESolver(derivative, "./include/initial-conditions.json", g_NSYS);

}
