

#include "decs.h"

/*

   given energy of photon in fluid rest frame w, in units of electron rest mass
   energy, and temperature of plasma, again in electron rest-mass units, return
   hot
   cross section in cgs.

   This has been checked against Wienke's Table 1, with some disagreement at
   the one part in 10^{-3} level, see wienke_table_1 in the subdirectory
   hotcross.
   It is not clear what this is due to, but Table 1 does appear to have been
   evaluated
   using Monte Carlo integration (!).

   A better way to do this would be to make a table in w*thetae and w/thetae;
   most
   of the variation is accounted for by w*thetae.

 */

#define MINW 1.e-12
#define MAXW 1.e9
#define MINT 1.e-3
#if KAPPA
/* A numerical table boundary, not a change to the physical cutoff. Valid
 * temperatures above this boundary use direct quadrature. */
#define MAXT fmin(TMAX, 0.5 * KAPPA_GAMMA_CUT)
#else
#define MAXT TMAX
#endif
#define NW 210
#define NT 100

#define KMONTY_STRINGIFY_IMPL(x) #x
#define KMONTY_STRINGIFY(x) KMONTY_STRINGIFY_IMPL(x)
#if KAPPA
#define HOTCROSS "hotcross_kappa_v2_" KMONTY_STRINGIFY(kappa_synch) "_" \
                 KMONTY_STRINGIFY(KAPPA_GAMMA_CUT) ".dat"
#else
#define HOTCROSS "hotcross.dat"
#endif

double table[NW + 1][NT + 1];
static double dlw, dlT, lminw, lmint;

static void hotcross_fail(const char *message) {
    fprintf(stderr, "kmonty hotcross: %s\n", message);
#if MPI
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
    exit(EXIT_FAILURE);
#endif
}

void init_hotcross(void) {
    int i, j, idum, jdum, nread;
    double lw, lT;
    double total_compton_cross_num(double w, double thetae);
    FILE *fp = NULL;
    int world_rank = 0;

    dlw = log10(MAXW / MINW) / NW;
    dlT = log10(MAXT / MINT) / NT;
    lminw = log10(MINW);
    lmint = log10(MINT);

#if MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
    /* Only rank zero touches the cache.  The table is broadcast below so a
     * multi-rank launch cannot race while creating or truncating the file. */
    if (world_rank == 0) {
        fp = fopen(HOTCROSS, "r");
        if (fp == NULL) {
            /* Reuse the matched width and normalization for each temperature. */
#pragma omp parallel for private(i, lw, lT) schedule(dynamic)
            for (j = 0; j <= NT; j++)
                for (i = 0; i <= NW; i++) {
                    double cross;
                    lw = lminw + i * dlw;
                    lT = lmint + j * dlT;
                    cross = total_compton_cross_num(pow(10., lw), pow(10., lT));
                    if (!(cross > 0.) || !isfinite(cross))
                        hotcross_fail("invalid generated cross section");
                    table[i][j] = log10(cross);
                }
            fp = fopen(HOTCROSS, "w");
            if (fp == NULL) hotcross_fail("could not write cache");
            for (i = 0; i <= NW; i++)
                for (j = 0; j <= NT; j++) {
                    lw = lminw + i * dlw;
                    lT = lmint + j * dlT;
                    if (fprintf(fp, "%d %d %.17g %.17g %.17g\n", i, j, lw, lT,
                                table[i][j]) < 0) hotcross_fail("cache write failed");
                }
        } else {
            for (i = 0; i <= NW; i++)
                for (j = 0; j <= NT; j++) {
                    nread = fscanf(fp, "%d %d %lf %lf %lf\n", &idum, &jdum,
                                   &lw, &lT, &table[i][j]);
                    if (nread != 5 || !isfinite(table[i][j]) ||
                        idum != i || jdum != j || !isfinite(lw) || !isfinite(lT) ||
                        fabs(lw-(lminw+i*dlw)) > 1.e-5 ||
                        fabs(lT-(lmint+j*dlT)) > 1.e-5)
                        hotcross_fail("invalid or truncated cache");
                }
        }
        if (fclose(fp) != 0) hotcross_fail("cache close failed");
    }

#if MPI
    MPI_Bcast(&table[0][0], (NW + 1) * (NT + 1), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
#endif
}

double total_compton_cross_lkup(double w, double thetae) {
    int i, j;
    double lw, lT, di, dj, lcross;
    double total_compton_cross_num(double w, double thetae);
    double hc_klein_nishina(double we);

#if KAPPA
    if (!(thetae > 0.) || !(thetae < KAPPA_GAMMA_CUT))
        hotcross_fail("temperature outside kappa matching domain");
#endif

    /* cold/low-energy: just use thomson cross section */
    if (w * thetae < 1.e-6)
        return (SIGMA_THOMSON);

    /* cold, but possible high energy photon: use klein-nishina */
    if (thetae < MINT)
        return (hc_klein_nishina(w) * SIGMA_THOMSON);

    /* in-bounds for table */
    if ((w > MINW && w < MAXW) && (thetae > MINT && thetae < MAXT)) {

        lw = log10(w);
        lT = log10(thetae);
        i = (int)((lw - lminw) / dlw);
        j = (int)((lT - lmint) / dlT);
        di = (lw - lminw) / dlw - i;
        dj = (lT - lmint) / dlT - j;

        lcross = (1. - di) * (1. - dj) * table[i][j] +
                 di * (1. - dj) * table[i + 1][j] +
                 (1. - di) * dj * table[i][j + 1] +
                 di * dj * table[i + 1][j + 1];

        if (isnan(lcross)) {
            fprintf(stderr, "lcross isnan %g %g %d %d %g %g\n", lw, lT, i, j,
                    di, dj);
        }

        return (pow(10., lcross));
    }

    fprintf(stderr, "out of bounds: %g %g\n", w, thetae);
    return (total_compton_cross_num(w, thetae));
}

#define MAXGAMMA 30.
#define DMUE 0.05
#define DGAMMAE 0.05

double calc_norm(double thetae) {
    // we integrate over beta since this is nicely bound betwoon 0 and 1,
    // indsted of 1 and infinity.
    double kappa_function_int(double beta, void *params);

    double norm;
    double result = 0., error;
    gsl_function integrand;
    integrand.function = &kappa_function_int;
    integrand.params = &thetae;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    /* Integrator parameters */
    double absolute_error = 0.;
    double relative_error = 1.e-6;
    size_t limit = 1000;

    double xmax = log(1. + 1000. * thetae);
    double width = kappa_width_matched(thetae);
    /* At low Thetae the κ core has width O(width) in ln(gamma).  A single
       adaptive call over the full tail can miss that core and return zero,
       making the subsequent normalized cross section evaluate 0*infinity. */
    double core = fmin(0.5, log1p(1024. * fmax(width, 1.e-12)));
    gsl_integration_qag(&integrand, 0., core, absolute_error, relative_error,
                        limit, GSL_INTEG_GAUSS61, w, &result, &error);
    if (core < xmax) {
        double tail = 0.;
        gsl_integration_qag(&integrand, core, xmax, absolute_error,
                            relative_error, limit, GSL_INTEG_GAUSS61, w,
                            &tail, &error);
        result += tail;
    }

    gsl_integration_workspace_free(w);
    norm = 1. / result;

    return norm;
}

double total_compton_cross_num(double w, double thetae) {
    double dmue, dgammae, mue, gammae, f, cross, norm;
    double dNdgammae(double thetae, double gammae);
    double boostcross(double w, double mue, double gammae);
    double hc_klein_nishina(double we);

    if (isnan(w)) {
        fprintf(stderr, "compton cross isnan: %g %g\n", w, thetae);
        return (0.);
    }

    /* check for easy-to-do limits */
    if (thetae < MINT && w < MINW)
        return (SIGMA_THOMSON);
    if (thetae < MINT)
        return (hc_klein_nishina(w) * SIGMA_THOMSON);

    dmue = DMUE;
    dgammae = thetae * DGAMMAE;

/* integrate over mu_e, gamma_e, where mu_e is the cosine of the
   angle between k and u_e, and the angle k is assumed to lie,
   wlog, along the z axis */
#if (KAPPA)
    norm = calc_norm(thetae);
#elif (THERMAL)
    norm = 1;
#endif
    cross = 0.;
    for (mue = -1. + 0.5 * dmue; mue < 1.; mue += dmue)
        for (gammae = 1. + 0.5 * dgammae; gammae < 1. + MAXGAMMA * thetae;
             gammae += dgammae) {

            f = 0.5 * norm * dNdgammae(thetae, gammae);

            cross += dmue * dgammae * boostcross(w, mue, gammae) * f;

            if (isnan(cross)) {
                fprintf(stderr, "cross is nan %g %g %g %g %g %g\n", w, thetae,
                        mue, gammae, dNdgammae(thetae, gammae),
                        boostcross(w, mue, gammae));
            }
        }

    return (cross * SIGMA_THOMSON);
}

double dNdgammae(double thetae, double gammae) {
    double dNdgammae_nth(double thetae, double gammae);
    double dNdgammae_th(double thetae, double gammae);

#if (KAPPA)
    return dNdgammae_nth(thetae, gammae);
#elif (THERMAL)
    return dNdgammae_th(thetae, gammae);
#endif
    return 0;
}

double dNdgammae_th(double thetae, double gammae) {
    double K2f;

    if (thetae > 1.e-2) {
        K2f = gsl_sf_bessel_Kn(2, 1. / thetae) * exp(1. / thetae);
    } else {
        K2f = sqrt(M_PI * thetae / 2.);
    }

    return ((gammae * sqrt(gammae * gammae - 1.) / (thetae * K2f)) *
            exp(-(gammae - 1.) / thetae));
}

// use beta instead of gamma, because its nicely bound between 0 and 1.
double kappa_function_int(double beta, void *params) {
    double gamma = exp(beta);

    double thetae = *(double *)params;
    double kappa = kappa_synch;
    double w = kappa_width_matched(thetae);
    if (!(w > 0.) || !isfinite(w)) return 0.;

    // double gamma = 1/sqrt(1-beta*beta);
    // return beta*beta* pow(gamma,5.)*pow(( 1 + ( gamma -
    // 1)/(kappa*thetae)),-(kappa+1))*exp(-gamma/gamma_max);

    return gamma * gamma * sqrt(gamma * gamma - 1.) *
           pow((1 + (gamma - 1) / (kappa * w)), -(kappa + 1)) *
           exp(-gamma / gamma_max);
}

double dNdgammae_nth(double thetae, double gammae) {
    double kappa = kappa_synch;
    double w = kappa_width_matched(thetae);
    if (!(w > 0.) || !isfinite(w)) return 0.;
    // This is not yet normalized!
    double dNdgam = gammae * sqrt(gammae * gammae - 1.) *
                    pow((1 + (gammae - 1) / (kappa * w)), -(kappa + 1)) *
                    exp(-gammae / gamma_max);
    return dNdgam;
}

double boostcross(double w, double mue, double gammae) {
    double we, boostcross, v;
    double hc_klein_nishina(double we);

    /* energy in electron rest frame */
    v = sqrt(gammae * gammae - 1.) / gammae;
    we = w * gammae * (1. - mue * v);

    boostcross = hc_klein_nishina(we) * (1. - mue * v);

    if (boostcross > 2) {
        fprintf(stderr, "w,mue,gammae: %g %g %g\n", w, mue, gammae);
        fprintf(stderr, "v,we, boostcross: %g %g %g\n", v, we, boostcross);
        fprintf(stderr, "kn: %g %g %g\n", v, we, boostcross);
    }

    if (isnan(boostcross)) {
        fprintf(stderr, "isnan: %g %g %g\n", w, mue, gammae);
        exit(0);
    }

    return (boostcross);
}

double hc_klein_nishina(double we) {
    double sigma;

    if (we < 1.e-3)
        return (1. - 2. * we);

    sigma = (3. / 4.) *
            (2. / (we * we) +
             (1. / (2. * we) - (1. + we) / (we * we * we)) * log(1. + 2. * we) +
             (1. + we) / ((1. + 2. * we) * (1. + 2. * we)));

    return (sigma);
}
