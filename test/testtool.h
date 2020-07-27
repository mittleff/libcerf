typedef struct {
    int failed;
    int total;
} result_t;

#define TEST_NEAR(result, limit, function_val, expected_val) \
    ++result.total; \
    if (fabs(function_val-expected_val)>limit*expected_val) { \
        printf("failure in subtest %i: %s\n", result.total, #function_val); \
        printf("- fct value %20.15g\n", function_val); \
        printf("- expected  %20.15g\n", expected_val); \
        printf("=> error %6.2g above limit %6.2g\n", \
               fabs(function_val-expected_val), limit*expected_val); \
        ++result.failed; \
    }
