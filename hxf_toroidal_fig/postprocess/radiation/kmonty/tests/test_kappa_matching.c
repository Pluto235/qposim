#include "../decs.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void) {
    double w1 = kappa_width_matched(0.1);
    double w2 = kappa_width_matched(1.0);
    assert(isfinite(w1) && w1 > 0.0);
    assert(isfinite(w2) && w2 > 0.0);
    assert(w1 != w2);
    puts("kappa energy matching: PASS");
    return 0;
}
