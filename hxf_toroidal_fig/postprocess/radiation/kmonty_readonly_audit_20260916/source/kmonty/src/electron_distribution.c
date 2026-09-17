#include "decs.h"
#include <float.h>
#include <string.h>

/* kappa_exp_v1: gamma_cut is a physical exponential scale, never an integration cap. */
#define KAPPA_MODEL_ID "kappa_exp_v1"
#define KAPPA_GAMMA_CUT 1.0e5
static _Thread_local double cached_theta = NAN, cached_kappa = NAN, cached_cut = NAN, cached_w = NAN;

static double logq(double g, double w, int kap) {
    if (!(g >= 1.0) || !(w > 0.0)) return -INFINITY;
    double l = log(g) + 0.5*log(fmax(g*g-1.0, DBL_MIN));
    if (kap) l += -(kappa_synch+1.0)*log1p((g-1.0)/(kappa_synch*w)) - g/KAPPA_GAMMA_CUT;
    else l += -g/w;
    return l;
}
static double moment(double theta, double w, int kap, int energy) {
    if (!(theta>0) || !(w>0)) return NAN;
    double xmax=log(2.0), prev=NAN;
    for (int pass=0; pass<16; ++pass) {
        int n=2048; double h=xmax/n, m=-INFINITY;
        for(int i=0;i<=n;i++){double x=i*h,g=exp(x),z=logq(g,w,kap);if(z>m)m=z;}
        double sum=0;
        for(int i=0;i<=n;i++){double x=i*h,g=exp(x),z=logq(g,w,kap);double v=exp(z-m)*(energy?(g-1):1)*g;sum+=(i==0||i==n?1:(i&1?4:2))*v;}
        double val=exp(m)*sum*h/3.;
        if (isfinite(prev) && fabs(val-prev)<=1e-11*fmax(1.,fabs(val))) return val;
        prev=val; xmax*=2.;
    }
    return prev;
}
static double mean(double theta,double w,int kap){double n=moment(theta,w,kap,0),e=moment(theta,w,kap,1);return n>0?e/n:NAN;}

double kappa_width_matched(double thetae) {
#if KAPPA
    if (!(thetae>0) || !(kappa_synch>1.0) || !isfinite(thetae)) return NAN;
    if (thetae==cached_theta && kappa_synch==cached_kappa && cached_cut==KAPPA_GAMMA_CUT) return cached_w;
    double target=mean(thetae,thetae,0); if(!isfinite(target)) return NAN;
    double lo=DBL_MIN*1e6, hi=fmax(1e-6,thetae), eh=mean(thetae,hi,1);
    for(int i=0;i<80 && (!isfinite(eh)||eh<target);i++){hi*=4.;eh=mean(thetae,hi,1);}
    if(!isfinite(eh)||eh<target) return NAN;
    for(int i=0;i<100;i++){double mid=sqrt(lo*hi), em=mean(thetae,mid,1);if(!isfinite(em)) return NAN;if(em<target)lo=mid;else hi=mid;}
    double w=sqrt(lo*hi); if(!isfinite(w)) return NAN;
    cached_theta=thetae;cached_kappa=kappa_synch;cached_cut=KAPPA_GAMMA_CUT;cached_w=w; return w;
#else
    (void)thetae; return NAN;
#endif
}
