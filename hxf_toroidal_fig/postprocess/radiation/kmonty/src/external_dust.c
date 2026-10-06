#include "external_dust.h"

extern double L_unit, T_unit, M_unit;
extern long int Ns;
extern gsl_rng *r;

static double dust_uniform(void) { return monty_rand(); }

void validate_electron_model(void) {
#if KAPPA
    if (!(kappa_synch > KAPPA_MIN_VALUE) ||
        !(KAPPA_GAMMA_CUT > gamma_min) ||
        !(gamma_max > gamma_min) || !isfinite(kappa_synch) ||
        !isfinite(gamma_max)) {
        fprintf(stderr, "invalid kappa configuration: kappa=%g gamma_min=%g gamma_max=%g TMAX=%g\n",
                kappa_synch, gamma_min, gamma_max, TMAX);
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
#endif
#if EXTERNAL_DUST
    {
        struct dust_field_config cfg;
        if (dust_field_configure(&cfg) != 1 ||
            !isfinite(dust_surface_flux(&cfg))) {
            fprintf(stderr, "invalid external dust configuration\n");
            MPI_Abort(MPI_COMM_WORLD, 3);
        }
    }
#endif
}

int dust_field_configure(struct dust_field_config *cfg) {
#if EXTERNAL_DUST
    if (!cfg) return -1;
    cfg->temperature_K = DUST_TEMPERATURE;
    cfg->luminosity_erg_s = DUST_COVERING * DUST_LDISK;
    cfg->radius_cm = DUST_RADIUS;
    cfg->center_offset_cm = DUST_CENTER_OFFSET;
    cfg->cross_section_cm = DUST_CROSS_SECTION * DUST_RADIUS;
    if (!(cfg->temperature_K > 0.0 && cfg->luminosity_erg_s > 0.0 &&
          cfg->radius_cm > 0.0 && cfg->cross_section_cm > 0.0)) return -1;
    return 1;
#else
    (void)cfg;
    return 0;
#endif
}

double dust_ring_area(const struct dust_field_config *cfg) {
    if (!cfg || cfg->radius_cm <= 0.0 || cfg->cross_section_cm <= 0.0)
        return 0.0;
    /* Surface area of the torus: 4 pi^2 R a. */
    return 4.0 * M_PI * M_PI * cfg->radius_cm * cfg->cross_section_cm;
}

double dust_surface_flux(const struct dust_field_config *cfg) {
    double area = dust_ring_area(cfg);
    return area > 0.0 ? cfg->luminosity_erg_s / area : 0.0;
}

const char *kmonty_seed_source_name(int source) {
    switch (source) {
    case KMONTY_SEED_SYNCHROTRON: return "synchrotron";
    case KMONTY_SEED_DUST_IR: return "dust_ir";
    default: return "unknown";
    }
}

static int dust_planck_frequency(double *nu_hz) {
    const double xmax = 40.0;
    const double ymax = 1.4215;
    for (int n = 0; n < 100000; ++n) {
        double x = xmax * dust_uniform();
        double y = ymax * dust_uniform();
        double f = x > 1.e-8 ? x*x*x / expm1(x) : 0.;
        if (y <= f) {
            *nu_hz = x * KBOL * DUST_TEMPERATURE / HPL;
            return isfinite(*nu_hz) && *nu_hz > 0. ? 0 : -1;
        }
    }
    return -1;
}

int dust_sample_packet(const double X[NDIM], double K[NDIM], double *nu_hz) {
#if EXTERNAL_DUST
    struct dust_field_config cfg;
    if (dust_field_configure(&cfg) != 1 || !X || !K || !nu_hz) return -1;
    double phi = 2. * M_PI * dust_uniform();
    double v = 2. * M_PI * dust_uniform();
    double a = cfg.cross_section_cm;
    double R = cfg.radius_cm;
    double source[3];
    source[0] = (R + a*cos(v))*cos(phi);
    source[1] = (R + a*cos(v))*sin(phi);
    source[2] = cfg.center_offset_cm + a*sin(v);
    double target[3] = {X[1]*L_unit, X[2]*L_unit, X[3]*L_unit};
    double d[3], norm2 = 0.;
    for (int q=0; q<3; ++q) { d[q] = target[q] - source[q]; norm2 += d[q]*d[q]; }
    if (!(norm2 > 0.) || !isfinite(norm2)) return -1;
    double inv = 1. / sqrt(norm2);
    if (dust_planck_frequency(nu_hz) != 0) return -1;
    double E = (*nu_hz) * HPL / (ME*CL*CL);
    K[0] = E;
    for (int q=0; q<3; ++q) K[q+1] = E*d[q]*inv;
    for (int q=0; q<4; ++q) if (!isfinite(K[q])) return -1;
    return 0;
#else
    (void)X; (void)K; (void)nu_hz;
    return -1;
#endif
}

double dust_packet_weight(double nu_hz, long int local_ns) {
#if EXTERNAL_DUST
    int world_size = 1;
#if MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
    if (!(nu_hz > 0.) || local_ns <= 0) return NAN;
    /* One local packet stream represents L_DT*T_unit/world_size.  Dividing by
     * h nu converts its energy packet into a photon-number weight. */
    double total = DUST_COVERING*DUST_LDISK*T_unit /
                   (ME*CL*CL) / (double)(local_ns*world_size*DUST_PACKET_FRACTION);
    return total / (nu_hz*HPL/(ME*CL*CL));
#else
    (void)nu_hz; (void)local_ns;
    return NAN;
#endif
}
