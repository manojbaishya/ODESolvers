/*
* @Author: manoj
* @Date:   2020-03-23 12:49:12
* @Last Modified by:   Manoj Baishya
* @Last Modified time: 2020-04-28 14:01:22
*/

#include "ODESolvers.h"
#include "derivatives.h"

const char gConfig[] = "./workspace/initial-conditions.json";

void singleODE(void);
void systemODE(void);

int main(int argc, char const *argv[]){

    singleODE();

    // systemODE();

    return 0;
}

void singleODE(void){

    set_parameters(&g_consts);
    callODESolver(derivative, events, gConfig, g_NSYS);

}

void systemODE(void){

    set_parameters(&g_consts);
    callODESolver(derivative, events, gConfig, g_NSYS);

}
