#ifndef KMONTY_EXTERNAL_DUST_H
#define KMONTY_EXTERNAL_DUST_H

#include "decs.h"

struct dust_field_config {
    double temperature_K;
    double luminosity_erg_s;
    double radius_cm;
    double center_offset_cm;
    double cross_section_cm;
};

/* Returns zero when the feature is disabled.  The implementation is a
 * geometry/normalisation contract only until the injection hook is enabled;
 * it never fabricates photons. */
int dust_field_configure(struct dust_field_config *cfg);
double dust_ring_area(const struct dust_field_config *cfg);
double dust_surface_flux(const struct dust_field_config *cfg);
const char *kmonty_seed_source_name(int source);

/* Sample one blackbody photon emitted by the finite torus and aimed at the
 * requested cell.  K is a contravariant AGN-frame null momentum in code
 * energy units. */
int dust_sample_packet(const double X[NDIM], double K[NDIM], double *nu_hz);
double dust_packet_weight(double nu_hz, long int local_ns);

#endif
