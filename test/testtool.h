typedef struct {
    int failed;
    int total;
} result_t;

#define TEST_NEAR(result, limit, function_val, expected_val)                                       \
    ++result.total;                                                                                \
    if (fabs(function_val - expected_val) > limit * fabs(expected_val)) {                          \
        printf("failure in subtest %i: %s\n", result.total, #function_val);                        \
        printf("- fct value %20.15g\n", function_val);                                             \
        printf("- expected  %20.15g\n", expected_val);                                             \
        printf(                                                                                    \
            "=> error %6.2g above limit %6.2g\n", fabs(function_val - expected_val),               \
            limit* fabs(expected_val));                                                            \
        ++result.failed;                                                                           \
    }

#define ZTEST(result, limit, function_val, expected_val)                                           \
    ++result.total;                                                                                \
    if (fabs(creal(function_val - expected_val)) > limit * fabs(creal(expected_val))               \
        || fabs(cimag(function_val - expected_val)) > limit * fabs(cimag(expected_val))) {         \
        printf("failure in subtest %i: %s\n", result.total, #function_val);                        \
        printf("- fct value %20.15g%+20.15gi\n", creal(function_val), cimag(function_val));        \
        printf("- expected  %20.15g%+20.15gi\n", creal(expected_val), cimag(expected_val));        \
        printf(                                                                                    \
            "=> error %6.2g or %6.2g above limit %6.2g or %6.2g \n",                               \
            fabs(creal(function_val - expected_val)), fabs(cimag(function_val - expected_val)),    \
            limit* fabs(creal(expected_val)), limit* fabs(cimag(expected_val)));                   \
        ++result.failed;                                                                           \
    }
