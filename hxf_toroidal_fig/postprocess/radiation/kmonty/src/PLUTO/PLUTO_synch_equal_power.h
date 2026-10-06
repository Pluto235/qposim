/* Private implementation included by PLUTO_utils.c. Explicit opt-in only.
 * Creation power in static Cartesian kappa4 cells; no transport changes. */
#include <errno.h>
#include <string.h>
#include <gsl/gsl_errno.h>

extern int kappa_F_log_grid(double, double, double *, double *);

struct synch_segment { double lo, hi, slope, cumulative; };
struct synch_cell {
    int n;
    double luminosity, integral, lambda;
    struct synch_segment segment[N_KAPPA_F];
};
static struct synch_cell *synch_cells;
static double Nsyn_global, synch_Lsum, synch_P0, synch_quad_error;

static size_t synch_index(int i, int j, int k) {
    return ((size_t)i * N2 + j) * N3 + k;
}

struct synch_integrand { double T, B; };
static double synch_power_integrand(double x, void *arg) {
    struct synch_integrand *p = arg;
    double nu = exp(x);
    return nu * F_eval_nth(p->T, p->B, nu);
}

static int init_synch_equal_power(void) {
    const char *mode = getenv("KMONTY_SYNCH_SOURCE");
    if (!mode || !strcmp(mode, "legacy")) return 0;
    if (strcmp(mode, "equal_power")) seed_error("unknown synch source", 0.);
    int ranks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (ranks != 1) {
        fprintf(stderr, "FATAL equal_power requires MPI=1\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (!KAPPA || THERMAL || MIXED || POWERLAW || kappa_synch != 4. ||
        COMPTON || EXTERNAL_DUST || !CAR || omp_get_max_threads() != 1)
        seed_error("equal_power requires static CAR kappa4, no Compton/dust, OMP=1", 0.);
    const char *budget = getenv("KMONTY_Nsyn_global");
    char *end;
    errno = 0;
    Nsyn_global = budget ? strtod(budget, &end) : 0.;
    if (!budget || end == budget || *end || errno || !isfinite(Nsyn_global) ||
        Nsyn_global <= 0. || Nsyn_global > INT_MAX-1)
        seed_error("invalid Nsyn_global", Nsyn_global);
    if (synch_cells) seed_error("equal_power initialized twice", 0.);
    synch_cells = calloc((size_t)N1*N2*N3, sizeof(*synch_cells));
    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);
    if (!synch_cells || !ws) seed_error("equal_power allocation", 0.);
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    for (int i=0; i<N1; i++) for (int j=0; j<N2; j++) for (int k=0; k<N3; k++) {
        double Ne, T, B, U[NDIM], Bc[NDIM], X[NDIM];
        int acc = 0;
        get_fluid_zone(i,j,k,&Ne,&T,&B,U,Bc,&acc);
        coord(i,j,k,X);
        if (hypot(X[1],X[2]) > RMAX || Ne == 0. || T < THETAE_MIN || B < 1e-20)
            continue;
        for (int a=0; a<NDIM; a++) {
            if (!isfinite(U[a]) || fabs(U[a]-(a==0)) > 1e-13)
                seed_error("equal_power requires static fluid", U[a]);
            for (int b=0; b<NDIM; b++)
                if (!isfinite(geom[i][j].gcov[a][b]) ||
                    fabs(geom[i][j].gcov[a][b]-(a==b?(a==0?-1.:1.):0.)) > 1e-13)
                    seed_error("equal_power requires Minkowski metric", geom[i][j].gcov[a][b]);
        }
        struct synch_cell *c = &synch_cells[synch_index(i,j,k)];
        double x[N_KAPPA_F+1], y[N_KAPPA_F+1];
        int intervals = kappa_F_log_grid(T,B,x,y);
        struct synch_integrand p = {T,B};
        gsl_function f = {synch_power_integrand, &p};
        double error = 0.;
        for (int q=0; q<intervals; q++) {
            double lo=fmax(log(NUMIN),x[q]), hi=fmin(log(NUMAX),x[q+1]);
            if (hi <= lo) continue;
            double value, err;
            /* Tiny cutoff segments hit floating-point roundoff in log(nu).
             * Allocate absolute error from accumulated power: at most
             * N_KAPPA_F*1e-14 of the final integral; retain the 1e-10
             * accumulated error gate below. No probability/weight change. */
            int status=gsl_integration_qag(&f,lo,hi,1e-14*c->integral,1e-11,1000,
                                           GSL_INTEG_GAUSS61,ws,&value,&err);
            /* The current angular kernel has its own hard cutoff. A -inf
             * log-table endpoint makes the open interpolation interval zero.
             * Preserve that existing support; never bridge it in the CDF. */
            if (!status && value==0. && err==0. &&
                (y[q]==-INFINITY || y[q+1]==-INFINITY)) continue;
            if (status || !isfinite(value) || value<=0.) {
                fprintf(stderr,"equal_power cell=%d,%d,%d segment=%d T=%.17g B=%.17g lognu=[%.17g,%.17g] integral=%.17g error=%.17g\n",
                        i,j,k,q,T,B,lo,hi,value,err);
                seed_error("equal_power segment quadrature", status);
            }
            c->integral += value;
            error += err;
            struct synch_segment *s = &c->segment[c->n++];
            s->lo=lo; s->hi=hi;
            s->slope=1.+(y[q+1]-y[q])/(x[q+1]-x[q]);
            s->cumulative=c->integral;
            if (!isfinite(s->slope)) seed_error("equal_power slope",s->slope);
        }
        double volume=dx[1][i]*dx[2][j]*dx[3][k]*pow(L_unit,3)*geom[i][j].g;
        c->luminosity=EE*EE*EE/(2.*M_PI*ME*CL*CL)*Ne*B*volume*c->integral;
        if (!isfinite(c->luminosity) || c->luminosity<0. || !isfinite(volume) || volume<=0.)
            seed_error("equal_power cell luminosity",c->luminosity);
        if (c->integral>0.) {
            synch_quad_error=fmax(synch_quad_error,error/c->integral);
            if (error/c->integral>1e-10) seed_error("equal_power quadrature error",error/c->integral);
        }
        synch_Lsum += c->luminosity;
    }
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(ws);
    if (!isfinite(synch_Lsum) || synch_Lsum<=0.) seed_error("empty equal_power source",synch_Lsum);
    synch_P0=synch_Lsum/Nsyn_global;
    for (size_t q=0; q<(size_t)N1*N2*N3; q++)
        synch_cells[q].lambda=Nsyn_global*(synch_cells[q].luminosity/synch_Lsum);
    fprintf(stderr,"equal_power Nsyn_global=%.17g Lsum=%.17g P0=%.17g quad_error=%.3g\n",
            Nsyn_global,synch_Lsum,synch_P0,synch_quad_error);
    return 1;
}

static double synch_equal_power_frequency(int i, int j, int k) {
    struct synch_cell *c=&synch_cells[synch_index(i,j,k)];
    if (!c->n || c->integral<=0.) seed_error("empty equal_power cell",c->integral);
    double target=monty_rand()*c->integral;
    int lo=0, hi=c->n-1;
    while (lo<hi) {
        int mid=(lo+hi)/2;
        if (target<c->segment[mid].cumulative) hi=mid; else lo=mid+1;
    }
    struct synch_segment *s=&c->segment[lo];
    /* In each real F interpolation interval, nu*F is exp(slope*x) times
     * exp(-nu/nu_cutoff). Invert the first factor and reject by the second.
     * Selecting intervals with their QAG masses gives the target density;
     * it does not silently substitute the uncorrected proposal. */
    for (int attempt=0; attempt<1000000; attempt++) {
        double u=monty_rand(), length=s->hi-s->lo, z=s->slope*length;
        double offset;
        if (fabs(z)<1e-10) offset=u*length;
        else if (z>50.) offset=length+log(u+(1.-u)*exp(-z))/s->slope;
        else offset=log1p(u*expm1(z))/s->slope;
        double x=fmax(nextafter(s->lo,s->hi),fmin(nextafter(s->hi,s->lo),s->lo+offset));
        double nu=exp(x);
        double p=probability_checked(exp(-(nu-exp(s->lo))/nu_cutoff));
        if (monty_rand()<p) return nu;
    }
    seed_error("equal_power frequency rejection exhausted",0.);
    return 0.;
}
