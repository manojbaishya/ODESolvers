#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector_long_double.h>
#include <gsl/gsl_matrix_long_double.h>


int main(int argc, char const *argv[]){

    gsl_matrix_long_double *mat = gsl_matrix_long_double_alloc(10, 10);

    long double num = 0;
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {

            num += 5;

            gsl_matrix_long_double_set(mat, i, j, pow(num, 2));

        }
    }

    gsl_vector_long_double_view viewCol = gsl_matrix_long_double_column(mat, 5);

    for (int i = 0; i < 10; ++i) {
        printf("Element %d is %Lf\n", i, gsl_vector_long_double_get(&viewCol.vector, i));
    }

    return 0;
}
