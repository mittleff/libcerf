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

IMPORT extern int faddeeva_algorithm;
IMPORT extern int faddeeva_nofterms;

int main( int argc, char **argv )
{
    double x, s, g;

    if( argc!=4 ){
        fprintf( stderr,  "usage:\n" );
        fprintf( stderr,  "   runvoigt <x> <sigma> <gamma>\n" );
        exit(-1);
    }

    x = atof( argv[1] );
    s = atof( argv[2] );
    g = atof( argv[3] );

    double y = voigt(x,s,g);
    printf( "%25.19g %3i %3i\n", y, faddeeva_algorithm, faddeeva_nofterms );
    return 0;
}
