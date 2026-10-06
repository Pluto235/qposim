#include "../decs.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void) {
    const double thetae[] = {0.01, 0.03, 0.1, 1.0, 3.0, 10.0, 100.0,
                             1.e3, 1.e4};
    const size_t nthetae = sizeof(thetae) / sizeof(thetae[0]);
    double first = 0.0;
    for (size_t i = 0; i < nthetae; ++i) {
        double w = kappa_width_matched(thetae[i]);
        assert(isfinite(w) && w > 0.0);
        if (i == 0) first = w;
        printf("kappa=%g Thetae=%g w=%.17g\n", kappa_synch, thetae[i], w);
    }
    assert(first != kappa_width_matched(thetae[nthetae - 1]));
    assert(isnan(kappa_width_matched(0.9 * KAPPA_GAMMA_CUT)));
    assert(isnan(kappa_width_matched(KAPPA_GAMMA_CUT)));
    assert(isnan(kappa_width_matched(2. * KAPPA_GAMMA_CUT)));
    puts("kappa energy matching grid: PASS");
    return 0;
}
