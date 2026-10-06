#include "../decs.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Deterministic probe; physics statistics are checked against the independent
 * CDF script, while this catches sampler failures and invalid support. */
double monty_rand(void) { return drand48(); }

int main(void) {
    const double temperatures[] = {0.01, 0.1, 1., 10., 100., 1.e3, 1.e4};
    srand48(1234);
    for (size_t t = 0; t < sizeof temperatures / sizeof temperatures[0]; ++t) {
        double sum = 0.;
        for (int i = 0; i < 1000; ++i) {
            double gamma, beta;
            sample_beta_distr_num(temperatures[t], &gamma, &beta);
            assert(isfinite(gamma) && gamma >= 1.);
            assert(isfinite(beta) && beta >= 0. && beta <= 1.);
            assert(gamma <= KAPPA_GAMMA_CUT * 64.);
            sum += gamma;
        }
        printf("kappa=%g Thetae=%g mean_gamma=%g PASS\n",
               kappa_synch, temperatures[t], sum / 1000.);
    }
    return 0;
}
