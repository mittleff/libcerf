/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File bigloop.c:
 *   Compute many function values, to measure timing.
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
#define size_t long int

int main()
{
    const size_t N = (size_t)1e9;
    const double step = 1./N;
    double sum = 0;

    for (size_t i=0; i<N; ++i)
        sum += im_w_of_x(i*step*40);
    printf("%g\n", sum/N);
    return 0;
}
