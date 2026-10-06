#ifndef DECS_H
#define DECS_H

#include "src/constants.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
//sphere model
//#define N1 (17288) // 512
//#define N2 (4096)  // 128
//#define N3 (1)
#ifndef COMPTON
#define COMPTON 1
#endif
//bhac & PLUTO model
extern int N1;
extern int N2;
extern int N3;

// Detector properies, amount of energy of angle bins.
#define N_ESAMP 130 //bins of super photons
#define N_EBINS 200 //bins of SED
#define N_THBINS 9
#define N_PHIBINS 1

#define rscreen 1e4 //radius of receiver screen

#define MPI 1
#define OPENMP 1

#define MKS 0
#define CKS 0
#define CAR 1

#define SFC 1
#define TRIM 0

#ifndef THERMAL
#define THERMAL 1
#endif
#ifndef KAPPA
#define KAPPA 0
#endif
#ifndef EXTERNAL_DUST
#define EXTERNAL_DUST 0
#endif
#define POWERLAW 0
#define MIXED 0

#ifndef kappa_synch
#define kappa_synch 4.0
#endif
#define perct_thermal 0.9
#define nu_cutoff 5.e33
#ifndef gamma_max
#define gamma_max 1.e5 //xf:carefully choose this value
#endif
#define gamma_min 1.
#define power_index 3.

/* External-field and electron-distribution contract.  These values are
 * deliberately compile-time for the first validated implementation: a run
 * cannot silently change its physics through an unrecorded environment. */
#define KAPPA_MIN_VALUE 3.000001
#define KAPPA_DEFAULT_VALUE 4.0
#define KAPPA8_VALUE 8.0
#define KAPPA_GAMMA_MAX 1.e5
#define DUST_TEMPERATURE 1000.0
#define DUST_LDISK 6.7e45
#define DUST_COVERING 0.6
#define DUST_RADIUS 6.5e18
#define DUST_CENTER_OFFSET 3.0e18
#define DUST_CROSS_SECTION 0.1
#define DUST_PACKET_FRACTION 0.5

#define NUMIN 1.e8
#define NUMAX 1.e20

#define FOLDING 0

#define DEBUG 0

#define NDIM 4
#define NPRIM 8

/* Range of initial superphoton frequencies */

#define THETAE_MIN 0.01
#define WEIGHT_MIN (1.e30)
#define TMIN (THETAE_MIN)
#define TMAX (1.e5)

#if (POWERLAW)
#define BTHSQMIN (1.e-4)
#define BTHSQMAX (1.e8)
#define NINT (10000)
#else
#define BTHSQMIN (1.e-4)
#define BTHSQMAX (1.e8)
#define NINT (10000)
#endif

/* mnemonics for primitive vars; conserved vars */
#define NP 8

#define KRHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7

#define D 0
#define S1 1
#define S2 2
#define S3 3
#define TAU 4
#define DS 8

/* numerical convenience */
#define SMALL 1.e-40

/* physical parameters */
#define MMW 0.5 /* mean molecular weight, in units of mp */

extern double MBH,R_LOW,R_HIGH;
extern double EPS,STEPMIN,STEPMAX,RMAX,SOURCE_DIST;
extern double trim_zi,trim_zf;
extern double weight_min;

/** data structures **/
struct of_photon {
    double X[NDIM];
    double K[NDIM];
    double dKdlam[NDIM];
    double w;
    double E;
    double L;
    double X1i;
    double X2i;
    double X3i;
    double tau_abs;
    double tau_scatt;
    double ne0;
    double thetae0;
    double b0;
    double E0;
    double E0s;
    int nscatt;
    /* 0: synchrotron, 1: external dust IR, future values reserved. */
    int seed_source;
};

enum kmonty_seed_source { KMONTY_SEED_SYNCHROTRON = 0,
                          KMONTY_SEED_DUST_IR = 1 };

struct of_geom {
    double gcon[NDIM][NDIM];
    double gcov[NDIM][NDIM];
    double g;
};

struct of_spectrum {
    double dNdlE;
    double dEdlE;
    double nph;
    double nscatt;
    double X1iav;
    double X2isq;
    double X3fsq;
    double tau_abs;
    double tau_scatt;
    double ne0;
    double thetae0;
    double b0;
    double E0;
    double source_scatt[2][3];
};

typedef struct of_spectrum type_spectr;

struct of_grid {
    struct of_spectrum spec[N_EBINS];
    double th, phi;
    int nlist;
    int *in;
};

struct block {
    int ind[3], level, size[3];
    double lb[3], dxc_block[3];
};

/** global variables **/
/** model independent */
extern gsl_rng *r;

extern double F_nth[N_ESAMP + 1], F_th[N_ESAMP + 1], wgt[N_ESAMP + 1];

extern long int Ns;
extern long int N_superph_recorded, N_superph_recorded_total, N_scatt;
extern int index;
/* HARM model globals */
extern struct of_geom **geom;
extern struct of_spectrum ***shared_spect;
extern double ***shared_Xi_spec;
extern double ***shared_ispec;
extern int n_within_horizon;

/* some coordinate parameters */
extern double a;
extern double R0, Rin, Rh, Rout, Rms;
extern double hslope;
extern double startx[NDIM], stopx[NDIM];
extern double dlE, lE0;
extern double gam;
extern double dMsim;

extern double M_unit;
extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Ne_unit;
extern double Thetae_unit;
extern double yhigh;

extern double max_tau_scatt, Ladv, dMact, bias_norm;

#define delta_num (1.e-7) // Used for numerical derivatives

/* some useful macros */
#define ZLOOP                                                                  \
    for (int i = 0; i < N1; i++)                                               \
        for (int j = 0; j < N2; j++)
#define DLOOP                                                                  \
    for (int k = 0; k < NDIM; k++)                                             \
        for (int l = 0; l < NDIM; l++)
#define MULOOP for (int mu = 0; mu < NDIM; mu++)
#define MUNULOOP                                                               \
    for (int mu = 0; mu < NDIM; mu++)                                          \
        for (int nu = 0; nu < NDIM; nu++)
#define LOOP_ijk                                                               \
    for (int i = 0; i < NDIM; i++)                                             \
        for (int j = 0; j < NDIM; j++)                                         \
            for (int k = 0; k < NDIM; k++)

#define INDEX(i, j, k) (NPRIM * ((k) + N3 * ((j) + N2 * (i))))

/** model-independent subroutines **/
/* core monte carlo/radiative transport routines */
void track_super_photon(struct of_photon *ph, long int *N_superph_recorded,
                        int igrid);
void record_super_photon(struct of_photon *ph);
void validate_electron_model(void);
double kappa_width_matched(double thetae);
void report_spectrum(double N_superph_made);
void scatter_super_photon(struct of_photon *ph, struct of_photon *php,
                          double Ne, double Thetae, double B, double Ucon[NDIM],
                          double Bcon[NDIM], double Gcov[NDIM][NDIM],
                          int ACCZONE);

/* OpenMP specific functions */
void omp_reduce_spect(void);

void mpi_reduce_spect(void);

/* MC/RT utilities */
void init_monty_rand(int seed);
double monty_rand(void);

/* geodesic integration */
void init_dKdlam(double X[], double Kcon[], double dK[]);
void push_photon_ham(double X[NDIM], double Kcon[][NDIM], double dl[]);
void push_photon(double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
                 double dl, double *E0, int n);
void push_photon4(double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
                  double dl);
void push_photon_cart(double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
                      double dl);
double stepsize(double X[NDIM], double K[NDIM], double dx_local);
void push_photon_gsl(double X[NDIM], double Kcon[NDIM], double dl);
int geodesic_deriv(double t, const double y[], double dy[], void *params);
void interpolate_geodesic(double Xi[], double X[], double Ki[], double K[],
                          double frac, double del_l);

/* basic coordinate functions supplied by grmonty */
void boost(double k[NDIM], double p[NDIM], double ke[NDIM]);
void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov);
double gdet_func(double gcov[][NDIM],
                 double X[NDIM]); /* calculated numerically */
void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
                          double K_tetrad[NDIM]);
void tetrad_to_coordinate(double Econ[NDIM][NDIM], double K_tetrad[NDIM],
                          double K[NDIM]);
double delta(int i, int j);
void normalize(double Ucon[NDIM], double Gcov[NDIM][NDIM]);
void normalize_null(double Gcov[NDIM][NDIM], double K[NDIM]);
void make_tetrad(double X[NDIM], double Ucon[NDIM], double Bhatcon[NDIM],
                 double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
                 double Ecov[NDIM][NDIM]);

/* functions related to basic radiation functions & physics */
/* physics-independent */
double get_fluid_nu(double X[4], double K[4], double Ucov[NDIM]);
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
                    double Bcov[NDIM], double B);
double alpha_inv_scatt(double nu, double thetae, double Ne);
double alpha_inv_abs(double nu, double thetae, double Ne, double B,
                     double theta, int ACCZONE );
double Bnu_inv(double nu, double thetae);
double jnu_inv(double nu, double thetae, double ne, double B, double theta);

double jnu_synch(double nu, double Ne, double Thetae, double B, double theta,
                 int ACCZONE);
double int_jnu(double Ne, double Thetae, double Bmag, double nu);
double jnu_integrand(double th, void *params);
double F_eval(double Thetae, double Bmag, double nu, int ACCZONE);
double jnu_bremss(double nu, double Ne, double Thetae);
double int_jnu_bremss(double Ne, double Thetae, double nu);

/* thermal synchrotron */
double jnu_synch_th(double nu, double Ne, double Thetae, double B,
                    double theta);
double int_jnu_th(double Ne, double Thetae, double Bmag, double nu);
void init_emiss_tables(void);
double F_eval_th(double Thetae, double Bmag, double nu);
double K2_eval(double Thetae);
double jnu_integrand_th(double th, void *params);
double alpha_inv_abs_th(double nu, double Thetae, double Ne, double B,
                        double theta);

/* non-thermal synchrotron */
double jnu_synch_nth(double nu, double Ne, double Thetae, double B,
                     double theta);
double int_jnu_nth(double Ne, double Thetae, double Bmag, double nu);
double anu_synch_nth(double nu, double Ne, double Thetae, double B,
                     double theta);
void init_emiss_tables_nth(void);
double jnu_integrand_nth(double th, void *params);
double F_eval_nth(double Thetae, double Bmag, double nu);
double alpha_inv_abs_nth(double nu, double Thetae, double Ne, double B,
                         double theta);
double get_kappa(double beta, double sigma);
double get_efficiency(double sigma, double beta);
double get_w_kappa(double theta_e, double beta, double sigma, double kappa);

/* compton scattering */
void init_hotcross(void);
double total_compton_cross_lkup(double nu, double theta);
double klein_nishina(double a, double ap);
double kappa_es(double nu, double theta);
void sample_electron_distr_p(double k[NDIM], double p[NDIM], double theta,
                             int ACCZONE );
void sample_beta_distr(double theta, double *gamma_e, double *beta_e,
                       int ACCZONE );
void sample_beta_distr_num(double theta, double *gamma_e, double *beta_e );
double sample_klein_nishina(double k0);
double sample_thomson(void);
double sample_mu_distr(double beta_e);
double sample_y_distr(double theta);
void sample_scattered_photon(double k[NDIM], double p[NDIM], double kp[NDIM]);

/*compton scattering kappa */
double sample_y_distr_nth(double Thetae);
double hypergeom_eval(double X);

/** model dependent functions required by code: these
   basic interfaces define the model **/

/* physics related */
void init_model(char *args[]);
void make_super_photon(struct of_photon *ph, int *quit_flag);
double bias_func(double Te, double w);
int get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
                     double *Thetae, double *B, double *sigma, double *beta,
                     double Ucon[NDIM], double Ucov[NDIM], double Bcon[NDIM],
                     double Bcov[NDIM], int *ACCZONE, double *dx_local,
                     int *igrid_c);

int stop_criterion(struct of_photon *ph);
int record_criterion(struct of_photon *ph);

/* coordinate related */

#endif
