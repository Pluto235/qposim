#include "../src/external_dust.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

double L_unit = 1.e15, T_unit = 1.e5, M_unit = 1.;
long int Ns = 1000;
gsl_rng *r = NULL;
double monty_rand(void) { return 0.1; }

int main(void) {
    struct dust_field_config c;
#if EXTERNAL_DUST
    assert(dust_field_configure(&c) == 1);
#else
    assert(dust_field_configure(&c) == 0);
#endif
    c.temperature_K = 1000.0;
    c.luminosity_erg_s = 0.6 * 6.7e45;
    c.radius_cm = 6.5e18;
    c.center_offset_cm = 3.0e18;
    c.cross_section_cm = 0.1 * c.radius_cm;
    assert(dust_ring_area(&c) > 0.0);
    assert(isfinite(dust_surface_flux(&c)) && dust_surface_flux(&c) > 0.0);
    assert(kmonty_seed_source_name(KMONTY_SEED_SYNCHROTRON)[0] == 's');
    assert(kmonty_seed_source_name(KMONTY_SEED_DUST_IR)[0] == 'd');
#if EXTERNAL_DUST
    {
        double X[NDIM] = {0., 0., 0., 0.};
        double K[NDIM], nu;
        assert(dust_sample_packet(X, K, &nu) == 0);
        assert(isfinite(nu) && nu > 0.0);
        for (int i=0; i<NDIM; ++i) assert(isfinite(K[i]));
    }
#endif
    puts("external dust contract: PASS");
    return 0;
}
