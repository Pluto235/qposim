#include "decs.h"
#include <float.h>

/* kappa_exp_v1: gamma_cut is a physical exponential scale, never an
 * integration cap.  The value is defined in decs.h so every κ module uses
 * the same cutoff. */
static _Thread_local double cached_theta = NAN, cached_kappa = NAN, cached_cut = NAN, cached_w = NAN;

static double logq(double g, double w, int kap) {
    if (!(g >= 1.0) || !(w > 0.0)) return -INFINITY;
    double l = log(g) + 0.5*log(fmax(g*g-1.0, DBL_MIN));
    if (kap) l += -(kappa_synch+1.0)*log1p((g-1.0)/(kappa_synch*w)) - g/KAPPA_GAMMA_CUT;
    else l += -(g-1.)/w;
    return l;
}
static double simpson_interval(double w, int kap, int energy, double a, double b, int n) {
    double h=(b-a)/n, m=-INFINITY;
    for(int i=0;i<=n;i++){
        double x=a+i*h,g=exp(x),z=logq(g,w,kap);
        if(z>m)m=z;
    }
    double sum=0;
    for(int i=0;i<=n;i++){
        double x=a+i*h,g=exp(x),z=logq(g,w,kap);
        double v=exp(z-m)*(energy?(g-1):1)*g;
        sum+=(i==0||i==n?1:(i&1?4:2))*v;
    }
    return exp(m)*sum*h/3.;
}
static double simpson_core(double w, int kap, int energy, double xmax, int n) {
    double h=sqrt(xmax)/n, m=-INFINITY;
    for(int i=0;i<=n;i++) {
        double z=i*h, x=z*z, g=exp(x), q=logq(g,w,kap);
        if(q>m) m=q;
    }
    double sum=0;
    for(int i=0;i<=n;i++) {
        double z=i*h, x=z*z, g=exp(x), q=logq(g,w,kap);
        double v=exp(q-m)*(energy?(g-1):1)*g*(2*z);
        sum+=(i==0||i==n?1:(i&1?4:2))*v;
    }
    return exp(m)*sum*h/3.;
}
static double thermal_moment(double theta, int energy) {
    const int n = 32768;
    const double zmax = 10.0;
    double h = zmax / n, sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        double z = i * h, y = z * z, g = 1.0 + theta * y;
        double v = 2.0 * z * theta * g * sqrt(fmax(g * g - 1.0, 0.0)) * exp(-y);
        if (energy) v *= theta * y;
        sum += (i == 0 || i == n ? 1.0 : (i & 1 ? 4.0 : 2.0)) * v;
    }
    return sum * h / 3.0;
}
static double moment(double theta, double w, int kap, int energy) {
    if (!(theta>0) || !(w>0)) return NAN;
#ifdef KMONTY_TEST_FORCE_QUADRATURE_FAILURE
    return NAN;
#endif
    if (!kap) return thermal_moment(theta, energy);
    /* The exponential cutoff makes contributions above this scale
       negligible, while the explicit bound keeps exp(x) finite during the
       adaptive refinement.  Without it, the old doubling loop eventually
       evaluated exp(x) at x beyond log(DBL_MAX) and returned NaN. */
    const double xmax_limit = log(KAPPA_GAMMA_CUT * 64.0);
    double xmax=log(2.0), prev=NAN;
    for (int pass=0; pass<16; ++pass) {
        /* Resolve the narrow non-relativistic core as well as the relativistic
           tail.  The previous 2048-point grid quantized low-Thetae kappa
           moments and could return the same width for distinct temperatures. */
        /* Resolve the non-relativistic core separately.  Its width is O(w),
           and a single uniform grid over the full tail misses it when w is
           much smaller than one. */
        double core=fmin(0.5, log1p(1024.*w));
        double val=simpson_core(w,kap,energy,core,8192);
        if (core < xmax) val += simpson_interval(w,kap,energy,core,xmax,16384);
        if (isfinite(prev) && fabs(val-prev)<=1e-11*fmax(1.,fabs(val))) return val;
        prev=val; xmax*=2.;
        if (xmax >= xmax_limit) {
            xmax = xmax_limit;
            if (pass == 15) break;
        }
    }
    return prev;
}
static double mean(double theta,double w,int kap){double n=moment(theta,w,kap,0),e=moment(theta,w,kap,1);return n>0?e/n:NAN;}

double kappa_width_matched(double thetae) {
#if KAPPA
    if (!(thetae>0) || !(kappa_synch>1.0) || !isfinite(thetae)) {
        fprintf(stderr, "kappa matching: INVALID_PARAMETER Thetae=%g kappa=%g\n", thetae, kappa_synch);
        return NAN;
    }
    if (thetae >= KAPPA_GAMMA_CUT) {
        fprintf(stderr, "kappa matching: NO_FINITE_ROOT Thetae=%g gamma_cut=%g\n", thetae, KAPPA_GAMMA_CUT);
        return NAN;
    }
    if (thetae==cached_theta && kappa_synch==cached_kappa && cached_cut==KAPPA_GAMMA_CUT) return cached_w;
    double target=mean(thetae,thetae,0);
    if(!isfinite(target)) {
        fprintf(stderr, "kappa matching: QUADRATURE_FAILED target Thetae=%g\n", thetae);
        return NAN;
    }
#ifdef KMONTY_TEST_FORCE_BRACKET_FAILURE
    fprintf(stderr, "kappa matching: BRACKET_FAILED injected Thetae=%g\n", thetae);
    return NAN;
#endif
    /* gamma_cut is the exponential scale, not a cap on the matched width. */
    double lo=1e-8, hi=fmax(1e-6,thetae);
    double eh=mean(thetae,hi,1);
    for(int i=0;i<80 && (!isfinite(eh)||eh<target) && hi<DBL_MAX/4.;i++){
        double next=hi*4.;
        if (next==hi) break;
        hi=next;
        eh=mean(thetae,hi,1);
    }
    if(!isfinite(eh)) {
        fprintf(stderr, "kappa matching: QUADRATURE_FAILED upper Thetae=%g\n", thetae);
        return NAN;
    }
    if(eh<target) {
        fprintf(stderr, "kappa matching: BRACKET_FAILED Thetae=%g target=%g upper=%g\n", thetae, target, eh);
        return NAN;
    }
    for(int i=0;i<100;i++){double mid=sqrt(lo*hi), em=mean(thetae,mid,1);if(!isfinite(em)) {fprintf(stderr, "kappa matching: QUADRATURE_FAILED bisection Thetae=%g\n", thetae); return NAN;}if(em<target)lo=mid;else hi=mid;}
    double w=sqrt(lo*hi); if(!isfinite(w) || !(w>0.)) return NAN;
    cached_theta=thetae;cached_kappa=kappa_synch;cached_cut=KAPPA_GAMMA_CUT;cached_w=w; return w;
#else
    (void)thetae; return NAN;
#endif
}

double kappa_width_required(double thetae) {
    double w = kappa_width_matched(thetae);
    if (!(w > 0.) || !isfinite(w)) {
        fprintf(stderr, "cannot match kappa electrons: Thetae=%g kappa=%g gamma_cut=%g\n",
                thetae, kappa_synch, KAPPA_GAMMA_CUT);
        exit(EXIT_FAILURE);
    }
    return w;
}
