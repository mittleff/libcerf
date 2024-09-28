/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File run_imwx.c:
 *   Interactive evaluation of Im w(x), which is proportional to dawson(x).
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
    fprintf(stderr, "   run_imwx x\n");
    exit(-1);
  }

  double x = atof(argv[1]);
  double y = im_w_of_x(x);

#ifdef CERF_INTROSPECT
  printf("%21.16e %22.17e %3i %3i\n", x, y, cerf_algorithm, cerf_nofterms);
#else
  printf("%21.16e %22.17e\n", x, y);
#endif
  return 0;
}
