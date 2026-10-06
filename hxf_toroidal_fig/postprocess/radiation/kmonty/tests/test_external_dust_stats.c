#include "../src/external_dust.h"
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>

double L_unit = 1.e15, T_unit = 1.e5, M_unit = 1.;
long int Ns = 1000;
gsl_rng *r = NULL;

static uint64_t state = UINT64_C(0x9e3779b97f4a7c15);
#pragma omp threadprivate(state)
double monty_rand(void) {
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return (double)((state * UINT64_C(2685821657736338717)) >> 11) /
           9007199254740992.0;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
#if EXTERNAL_DUST
    struct dust_field_config c;
    assert(dust_field_configure(&c) == 1);
    double area = 4. * M_PI * M_PI * c.radius_cm * c.cross_section_cm;
    assert(fabs(dust_ring_area(&c) / area - 1.) < 1.e-14);
    assert(fabs(dust_surface_flux(&c) * area / c.luminosity_erg_s - 1.) < 1.e-14);

    const double X[NDIM] = {0., 0., 0., 0.};
    const int n = 20000;
    double mean_x = 0., mean_ux = 0., mean_uy = 0., mean_ux2 = 0., mean_uy2 = 0.;
    #pragma omp parallel for reduction(+:mean_x,mean_ux,mean_uy,mean_ux2,mean_uy2)
    for (int i = 0; i < n; ++i) {
        state = UINT64_C(0x9e3779b97f4a7c15) ^
                (uint64_t)(i + 1) * UINT64_C(0x94d049bb133111eb);
        double K[NDIM], nu;
        assert(dust_sample_packet(X, K, &nu) == 0);
        assert(isfinite(nu) && nu > 0.);
        assert(K[0] > 0. && isfinite(K[0]));
        double ux = K[1] / K[0], uy = K[2] / K[0], uz = K[3] / K[0];
        assert(fabs(ux * ux + uy * uy + uz * uz - 1.) < 1.e-12);
        /* The target is below the torus, so every ray must point toward -z. */
        assert(uz < 0.);
        double x = nu * HPL / (KBOL * c.temperature_K);
        mean_x += x;
        mean_ux += ux;
        mean_uy += uy;
        mean_ux2 += ux * ux;
        mean_uy2 += uy * uy;
    }
    mean_x /= n;
    mean_ux /= n;
    mean_uy /= n;
    mean_ux2 /= n;
    mean_uy2 /= n;
    /* Planck number-spectrum mean x is  pi^4/(30 zeta(3)) = 3.832. */
    assert(mean_x > 3.75 && mean_x < 3.92);
    assert(fabs(mean_ux) < 0.02 && fabs(mean_uy) < 0.02);
    assert(fabs(mean_ux2 - mean_uy2) < 0.02);

    int world_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    double K[NDIM], nu;
    assert(dust_sample_packet(X, K, &nu) == 0);
    double weight = dust_packet_weight(nu, Ns);
    double expected_energy = c.luminosity_erg_s * T_unit / world_size;
    double packet_energy = weight * nu * HPL * Ns * DUST_PACKET_FRACTION;
    assert(fabs(packet_energy / expected_energy - 1.) < 1.e-12);
    printf("dust stats: mean_x=%.6f mean_ux=%.6g mean_uy=%.6g "
           "var_xy_delta=%.6g PASS\n", mean_x, mean_ux, mean_uy,
           mean_ux2 - mean_uy2);
#else
    assert(dust_field_configure(NULL) == 0);
    puts("dust stats: EXTERNAL_DUST=0 PASS");
#endif
    MPI_Finalize();
    return 0;
}
