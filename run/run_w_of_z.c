/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File runvoigt.c:
 *   Interactive evaluation of Voigt's function.
 *
 * Copyright:
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Public domain.
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 */

#include <stdio.h>
#include <stdlib.h>
#include "cerf.h"
#include "defs.h"

IMPORT extern int faddeeva_algorithm;
IMPORT extern int faddeeva_nofterms;

int main( int argc, char **argv )
{
    double x, y;

    if( argc!=3 ){
        fprintf( stderr,  "usage:\n" );
        fprintf( stderr,  "   run_w_of_z <x> <y>\n" );
        exit(-1);
    }

    x = atof( argv[1] );
    y = atof( argv[2] );

    _cerf_cmplx w = w_of_z( C(x,y) );

    double v[2][2];
    v[0][0] = creal(w);
    v[0][1] = cimag(w);

    printf( "%25.19g %25.19g %3i %3i\n", v[0][0], v[0][1],
            faddeeva_algorithm, faddeeva_nofterms );
    return 0;
}
