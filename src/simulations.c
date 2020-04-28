/*
* @Author: manoj
* @Date:   2020-03-23 12:49:12
* @Last Modified by:   Manoj Baishya
* @Last Modified time: 2020-04-28 14:01:22
*/
#include "simulations.h"
#include "ODESolvers.h"
#include "derivatives.h"

const char gConfig[] = "./include/initial-conditions.json";

void singleODE(void){

    set_parameters(&g_consts);
    callODESolver(derivative, events, gConfig, g_NSYS);

}

void systemODE(void){

    set_parameters(&g_consts);
    callODESolver(derivative, events, gConfig, g_NSYS);

}
