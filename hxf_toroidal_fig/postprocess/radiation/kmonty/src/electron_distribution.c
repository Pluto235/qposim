#include "decs.h"
#include <float.h>

#define MOMENT_BINS 1024

static double cached_theta = -1.;
static double cached_w = NAN;
#pragma omp threadprivate(cached_theta, cached_w)

static double distribution_mean_kinetic(double theta, double width,
                                        int use_kappa) {
    double logmax = log(gamma_max);
    double norm = 0., energy = 0.;
    for (int i=0; i<MOMENT_BINS; ++i) {
        double lg = (i + 0.5) * logmax / MOMENT_BINS;
        double g = exp(lg);
        double base = g*g*sqrt(fmax(g*g-1., 0.)); /* dgamma=g dlg */
        double shape;
        if (use_kappa)
            shape = pow(1. + (g-1.)/(kappa_synch*width),
                        -(kappa_synch+1.)) * exp(-g/gamma_max);
        else
            shape = exp(-(g-1.)/theta);
        double q = base*shape;
        norm += q;
        energy += q*(g-1.);
    }
    return norm > 0. && isfinite(norm) && isfinite(energy) ? energy/norm : NAN;
}

double kappa_width_matched(double thetae) {
#if KAPPA
    if (!(thetae > 0.) || !(kappa_synch > KAPPA_MIN_VALUE) ||
        !(gamma_max > 1.)) return NAN;
    if (thetae == cached_theta) return cached_w;
    double target = distribution_mean_kinetic(thetae, thetae, 0);
    if (!isfinite(target)) return NAN;
    double lo = fmax(DBL_MIN, thetae*1.e-8), hi = thetae*1.e4;
    double ehi = distribution_mean_kinetic(thetae, hi, 1);
    if (!isfinite(ehi) || ehi < target) return NAN;
    for (int n=0; n<80; ++n) {
        double mid = sqrt(lo*hi);
        double emid = distribution_mean_kinetic(thetae, mid, 1);
        if (!isfinite(emid)) return NAN;
        if (emid < target) lo = mid; else hi = mid;
    }
    double w = sqrt(lo*hi);
    double achieved = distribution_mean_kinetic(thetae, w, 1);
    if (!isfinite(w) || fabs(achieved-target) > 1.e-6*fmax(target,1.e-12))
        return NAN;
    cached_theta = thetae;
    cached_w = w;
    return w;
#else
    (void)thetae;
    return NAN;
#endif
}
