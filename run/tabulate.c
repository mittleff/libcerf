/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File tabulate.c:
 *   Tabulate outcomes. Also used to generate test cases.
 *
 * Copyright:
 *   (C) 2022 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Public domain.
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, 2022
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 */

#include "cerf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void tabulate(double x) {
    printf( "    RTEST(result, 1e-13, im_w_of_x(%25.15g), %25.15g);\n", x, im_w_of_x(x) );
}

int main()
{
    tabulate(0);
    printf("\n    // rough logarithmic grid\n");
    for (int i=-275; i<=275; i += 50) {
        double x = pow(10., i);
        tabulate(-x);
        tabulate(x);
    }
    printf("\n    // medium logarithmic grid\n");
    for (int i=-15; i<=15; i += 2) {
        double x = pow(10., i);
        tabulate(-x);
        tabulate(x);
    }
    printf("\n    // fine logarithmic grid\n");
    for (double i=-3; i<=3; i += 0.2) {
        double x = pow(10., i);
        tabulate(-x);
        tabulate(x);
    }
    return 0;
}
