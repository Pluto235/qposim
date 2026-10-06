#include "decs.h"
#include <float.h>

/*
 * Direct rejection sampler for the κ electron distribution.  The old
 * implementation found an envelope maximum with an unchecked Brent bracket;
 * at high temperature GSL could abort before a photon was produced. Bound
 * each log-gamma interval using its increasing phase space factor at the
 * right edge and decreasing tail factor at the left edge.
 */

static double log_fdist(double ge, double w) {
#if KAPPA
    double kap = kappa_synch;
    if (!(ge > 1.0) || !(w > 0.0) || !isfinite(ge) || !isfinite(w))
        return -INFINITY;
    return 2.0 * log(ge) +
           0.5 * log(fmax((ge - 1.0) * (ge + 1.0), DBL_MIN)) -
           (kap + 1.0) * log1p((ge - 1.0) / (kap * w)) -
           ge / KAPPA_GAMMA_CUT;
#else
    (void)ge;
    (void)w;
    return -INFINITY;
#endif
}

static _Thread_local double envelope_theta = NAN;
static _Thread_local double envelope_width = NAN;
static _Thread_local double envelope_hi = NAN;
static _Thread_local double envelope_log = -INFINITY;

/* dN/d log gamma, retained as a small public probe for diagnostics. */
double fdist(double ge, double Thetae) {
#if KAPPA
    double w = kappa_width_matched(Thetae);
    if (!(w > 0.0) || !isfinite(w)) return 0.0;
    double lf = log_fdist(ge, w);
    return isfinite(lf) ? exp(lf) : 0.0;
#else
    if (!(ge > 1.0) || !(Thetae > 0.0)) return 0.0;
    return ge * ge * sqrt(fmax(ge * ge - 1.0, 0.0)) *
           exp(-(ge - 1.0) / Thetae);
#endif
}

void sample_beta_distr_num(double Thetae, double *gamma_e, double *beta_e) {
#if KAPPA
    enum { N_GRID = 8192, MAX_TRIES = 10000000 };
    double w = kappa_width_matched(Thetae);
    if (!(w > 0.0) || !isfinite(w)) {
        fprintf(stderr, "invalid energy-matched kappa width: thetae=%g\n",
                Thetae);
        exit(EXIT_FAILURE);
    }

    /* The exp(-gamma/gamma_cut) tail beyond 64 cutoffs is negligible. */
    double ge_lo = 1.0;
    double ge_hi = KAPPA_GAMMA_CUT * 64.0;
    if (!(ge_hi > ge_lo) || !isfinite(ge_hi)) {
        fprintf(stderr, "invalid kappa sampler window: thetae=%g w=%g\n",
                Thetae, w);
        exit(EXIT_FAILURE);
    }

    if (Thetae != envelope_theta || w != envelope_width) {
        double lge_lo = log(ge_lo), lge_hi = log(ge_hi);
        double log_fmax = -INFINITY;
        for (int i = 0; i < N_GRID; ++i) {
            double a = exp(lge_lo + (double)i / N_GRID * (lge_hi - lge_lo));
            double b = exp(lge_lo + (double)(i+1) / N_GRID * (lge_hi - lge_lo));
            double lf = 2.*log(b) + 0.5*log((b-1.)*(b+1.)) -
                (kappa_synch+1.)*log1p((a-1.)/(kappa_synch*w)) - a/KAPPA_GAMMA_CUT;
            if (lf > log_fmax) log_fmax = lf;
        }
        if (!isfinite(log_fmax)) {
            fprintf(stderr, "invalid kappa sampler envelope: thetae=%g w=%g\n",
                    Thetae, w);
            exit(EXIT_FAILURE);
        }
        envelope_theta = Thetae;
        envelope_width = w;
        envelope_hi = ge_hi;
        envelope_log = log_fmax + 1.e-12;
#ifdef KMONTY_TEST_FORCE_ENVELOPE_VIOLATION
        envelope_log -= 100.0;
#endif
    }

    double lge_lo = log(ge_lo), lge_hi = log(envelope_hi);
    /* An envelope violation is a numerical failure, never a clipped draw. */
    for (int n = 0; n <
#ifdef KMONTY_TEST_FORCE_SAMPLER_EXHAUSTION
             0
#else
             MAX_TRIES
#endif
             ; ++n) {
        double u = monty_rand();
        if (u <= 0.0) u = DBL_MIN;
        double ge = exp(lge_lo + monty_rand() * (lge_hi - lge_lo));
        double log_accept = log_fdist(ge, w) - envelope_log;
        if (log_accept > 1.e-12) {
            fprintf(stderr, "kappa sampler envelope violation: thetae=%g w=%g\n",
                    Thetae, w);
            exit(EXIT_FAILURE);
        }
        if (log(u) <= log_accept) {
            *gamma_e = ge;
            *beta_e = sqrt(fmax(0.0, 1.0 - 1.0 / (ge * ge)));
            return;
        }
    }

    fprintf(stderr, "kappa electron rejection sampler exhausted: thetae=%g w=%g\n",
            Thetae, w);
    exit(EXIT_FAILURE);
#else
    (void)Thetae;
    (void)gamma_e;
    (void)beta_e;
    fprintf(stderr, "kappa electron sampler called without KAPPA=1\n");
    exit(EXIT_FAILURE);
#endif
}
