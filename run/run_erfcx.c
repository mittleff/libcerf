/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File run_erfcx.c:
 *   Interactive evaluation of erfcx(x).
 *
 * Copyright:
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Public domain.
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, 2024
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 */

#include "cerf.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef CERF_INTROSPECT
IMPORT extern int cerf_algorithm;
IMPORT extern int cerf_nofterms;
#endif

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage:\n");
    fprintf(stderr, "   run_erfcx x\n");
    exit(-1);
  }

  double x = atof(argv[1]);
  double y = erfcx(x);

#ifdef CERF_INTROSPECT
  printf("%21.16e %22.17e %3i %3i\n", x, y, cerf_algorithm, cerf_nofterms);
#else
  printf("%21.16e %22.17e\n", x, y);
#endif
  return 0;
}
