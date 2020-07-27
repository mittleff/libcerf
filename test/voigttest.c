#include "cerf.h"
#include "testtool.h"
#include <math.h>
#include <stdio.h>

int main() {
    result_t result = {0,0};

    // expected results analytically determined:
    TEST_NEAR(result, 1e-15, voigt(0,1,0), 1/sqrt(6.283185307179586));
    TEST_NEAR(result, 1e-15, voigt(0,0,1), 1/3.141592653589793);
    TEST_NEAR(result, 1e-13, voigt(0,.5,.5), .41741856104074);

    // expected results obtained from scipy.integrate:
    TEST_NEAR(result, 1e-12, voigt(1,.5,.5), .18143039885260323);
    TEST_NEAR(result, 1e-12, voigt(1e5,.5e5,.5e5), .18143039885260323e-5);
    TEST_NEAR(result, 1e-12, voigt(1e-5,.5e-5,.5e-5), .18143039885260323e5);
    TEST_NEAR(result, 1e-12, voigt(1,.2,5), 0.06113399719916219);
    TEST_NEAR(result, 1e-12, voigt(1,5,.2), 0.07582140674553575);

    printf("%i/%i tests failed", result.failed, result.total);
    return result.failed ? 1 : 0;
}
