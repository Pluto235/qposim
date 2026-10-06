/* Private direct continuum backend, included only by jnu_mixed.c.
 * Radiation identity: kappa_continuum_gsl_v1; electron identity: kappa_exp_v1.
 * N(gamma)=Ne*q/Z, q=gamma*sqrt(gamma^2-1)*[1+(gamma-1)/(k*w)]^(-k-1)
 * * exp(-gamma/gcut). Isotropic phase density N/[4pi m^3 c^3 gamma p].
 * Pnu=sqrt(3)e^3 B sin(theta)/(mc^2) S(nu/[1.5 nuB sin(theta) gamma^2]);
 * j=Ne/(4pi)<Pnu>, alpha=Ne/(8pi m nu^2)<Pnu*((k+1)/(kw+gamma-1)+1/gcut)>.
 * Continuum/beaming approximation only, NOT a harmonic or quantum solution.
 * No legacy K/X cutoffs. Photon exp(-nu/nu_cutoff) is a separate common factor.
 */
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_errno.h>
#ifndef KAPPA_SYNCH_ELECTRON_ORDER
#define KAPPA_SYNCH_ELECTRON_ORDER 512
#endif
#ifndef KAPPA_SYNCH_ANGLE_ORDER
#define KAPPA_SYNCH_ANGLE_ORDER 64
#endif
#define SYNCH_NUM_ID "kappa_exp_v1/kappa_continuum_gsl_v1"
#define SYNCH_F_INTERP_TOL 5.e-4
static void synch_num_fail(const char *where, double value) {
    fprintf(stderr,"FATAL %s %s value=%.17g\n", SYNCH_NUM_ID,where,value);
    int initialized=0; MPI_Initialized(&initialized);
    if(initialized) MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    exit(EXIT_FAILURE);
}
static void synch_num_serial(void) {
    int ranks=1,initialized=0;MPI_Initialized(&initialized);
    if(initialized) MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    if(ranks!=1 || omp_get_max_threads()!=1)
        synch_num_fail("direct backend currently requires MPI=1 OMP=1",ranks);
}
/* One-entry in-memory electron quadrature cache. Key T and matched width;
 * kappa, gamma_cut, quadrature order, upper multiplier and algorithm identity
 * are immutable compile-time properties of this executable. No disk cache. */
static double sn_T=-1., sn_width, sn_tail_probability_bound, sn_g[KAPPA_SYNCH_ELECTRON_ORDER];
static double sn_f[KAPPA_SYNCH_ELECTRON_ORDER],sn_d[KAPPA_SYNCH_ELECTRON_ORDER];
static void synch_num_electrons(double T) {
    synch_num_serial();
    if(!isfinite(T)||T<=0.)synch_num_fail("invalid temperature",T);
    double width=kappa_width_required(T);
    if(T==sn_T && width==sn_width)return;
    double upper=100.*KAPPA_GAMMA_CUT, end=sqrt(log(upper)), Z=0.;
    gsl_integration_glfixed_table *tab=gsl_integration_glfixed_table_alloc(KAPPA_SYNCH_ELECTRON_ORDER);
    if(!tab)synch_num_fail("electron allocation",0.);
    for(int i=0;i<KAPPA_SYNCH_ELECTRON_ORDER;i++) {
        double t,v;gsl_integration_glfixed_point(0.,end,i,&t,&v,tab);
        double g=exp(t*t),p=sqrt(expm1(t*t)*(g+1.));
        double q=g*p*exp(-(kappa_synch+1.)*log1p((g-1.)/(kappa_synch*width))-g/KAPPA_GAMMA_CUT);
        sn_g[i]=g;sn_f[i]=v*2.*t*g*q;Z+=sn_f[i];
        sn_d[i]=(kappa_synch+1.)/(kappa_synch*width+g-1.)+1./KAPPA_GAMMA_CUT;
    }
    gsl_integration_glfixed_table_free(tab);
    if(!isfinite(Z)||Z<=0.)synch_num_fail("normalization",Z);
    for(int i=0;i<KAPPA_SYNCH_ELECTRON_ORDER;i++)sn_f[i]/=Z;
    /* For gamma>=U, q<=C gamma^(1-k) exp(-gamma/gcut),
     * C=[kw/(1-1/U)]^(k+1); integrate the decreasing power by its
     * endpoint upper bound. S(x)<1 makes this a coefficient tail bound. */
    sn_tail_probability_bound=exp((kappa_synch+1.)*log(kappa_synch*width/(1.-1./upper))
        +(1.-kappa_synch)*log(upper)+log(KAPPA_GAMMA_CUT)-upper/KAPPA_GAMMA_CUT-log(Z));
    sn_T=T;sn_width=width;
}
static void synch_num_coeff(double nu,double Ne,double T,double B,double theta,double *j,double *a) {
    if(!isfinite(nu)||nu<=0.||!isfinite(Ne)||Ne<0.||!isfinite(B)||B<0.||
       !isfinite(theta)||theta<0.||theta>M_PI)synch_num_fail("invalid coefficient input",nu);
    synch_num_electrons(T);
    *j=*a=0.;
    if(Ne==0.||B==0.||theta==0.||theta==M_PI)return;
    double s=sin(theta),nub=EE*B/(2.*M_PI*ME*CL),js=0.,as=0.;
    for(int i=0;i<KAPPA_SYNCH_ELECTRON_ORDER;i++) {
        double x=nu/(1.5*nub*s*sn_g[i]*sn_g[i]);
        /* S(x)>0 here, but S(x>700)<3.3e-303. Since sn_f sums to
         * one, the omitted j is bounded by pref*3.3e-303/(4pi),
         * and alpha by the same bound times max sn_d/(2m nu^2). */
        if(x>700.)continue;
        gsl_sf_result kernel;
        int status=gsl_sf_synchrotron_1_e(x,&kernel);
        /* GSL explicitly reports floating underflow. This is not a physical
         * support boundary; the omitted kernel is below representable range. */
        if(status==GSL_EUNDRFLW)continue;
        if(status||!isfinite(kernel.val)||kernel.val<0.)synch_num_fail("synchrotron kernel",status);
        double v=sn_f[i]*kernel.val;js+=v;as+=v*sn_d[i];
    }
    double omitted=sn_tail_probability_bound+3.3e-303;
    double derivative_at_upper=(kappa_synch+1.)/(kappa_synch*sn_width+100.*KAPPA_GAMMA_CUT-1.)+1./KAPPA_GAMMA_CUT;
    /* Include the change of normalization as well as the missing tail. */
    double dmax=(kappa_synch+1.)/(kappa_synch*sn_width)+1./KAPPA_GAMMA_CUT;
    double omitted_abs=sn_tail_probability_bound*fmax(as,derivative_at_upper)+3.3e-303*dmax;
    if(!(js>0.)||!(as>0.)||omitted/js>1e-7||omitted_abs/as>1e-7) {
        /* Near the field axis a relative error is meaningless for vanishing
         * coefficients. Certify the absolute error against the same nu,T,B
         * transverse coefficient, never against a fitted or rescaled source.
         * At the transverse direction itself, an unresolved tail is fatal.
         * Thus a far-out-of-range frequency cannot quietly become zero. */
        if(theta==M_PI/2.)synch_num_fail("coefficient outside controlled numerical tail/range",nu);
        double j90,a90,pref90=sqrt(3.)*EE*EE*EE*B/(ME*CL*CL);
        synch_num_coeff(nu,Ne,T,B,M_PI/2.,&j90,&a90);
        double jscale=j90*4.*M_PI/(Ne*pref90);
        double ascale=a90*8.*M_PI*ME*nu*nu/(Ne*pref90);
        if(omitted*s/jscale>1e-7||omitted_abs*s/ascale>1e-7)
            synch_num_fail("angle absolute tail error",nu);
    }
    double pref=sqrt(3.)*EE*EE*EE*B*s/(ME*CL*CL);
    *j=Ne*pref*js/(4.*M_PI);
    *a=Ne*pref*as/(8.*M_PI*ME*nu*nu);
    if(!isfinite(*j)||!isfinite(*a))synch_num_fail("coefficient overflow",nu);
}
double kappa_synch_num_j(double nu,double Ne,double T,double B,double theta) {
    double j,a;synch_num_coeff(nu,Ne,T,B,theta,&j,&a);
    return j*exp(-nu/nu_cutoff);
}
double kappa_synch_num_alpha(double nu,double Ne,double T,double B,double theta) {
    double j,a;synch_num_coeff(nu,Ne,T,B,theta,&j,&a);
    return a; /* alpha_inv_abs applies the common photon factor once. */
}
/* F before photon factor. Jnu=4pi int_0^1 j(mu)dmu=JCST*Ne*B*F.
 * This calls precisely the same electron/kernel routine as public jnu. */
double kappa_synch_num_F_direct(double T,double B,double nu) {
    if(!isfinite(B)||B<=0.)synch_num_fail("F requires positive B",B);
    gsl_integration_glfixed_table *tab=gsl_integration_glfixed_table_alloc(KAPPA_SYNCH_ANGLE_ORDER);
    if(!tab)synch_num_fail("angle allocation",0.);
    double J=0.;
    for(int i=0;i<KAPPA_SYNCH_ANGLE_ORDER;i++) {
        double theta,v,j,a;gsl_integration_glfixed_point(0.,M_PI/2.,i,&theta,&v,tab);
        synch_num_coeff(nu,1.,T,B,theta,&j,&a);J+=v*j*sin(theta);
    }
    gsl_integration_glfixed_table_free(tab);
    return 4.*M_PI*J/(EE*EE*EE*B/(2.*M_PI*ME*CL*CL));
}
/* Optional one-cell frequency representation needed by the existing equal-
 * power source. NOT a universal K table. Rebuilt for T,B,width changes; all
 * other model/precision/identity fields are immutable in this executable.
 * Outside [NUMIN,NUMAX], F is evaluated directly, never extrapolated or zero.
 * Midpoint checks bound the tested interpolation errors, not an analytic
 * supremum proof. The independent validation also checks off-midpoint nodes.
 */
static int sn_grid_ready;
static double sn_grid_T,sn_grid_B,sn_grid_width;
static double sn_x[N_KAPPA_F+1],sn_y[N_KAPPA_F+1];
static int synch_num_grid(double T,double B,double *x,double *y) {
    synch_num_electrons(T);
    if(!sn_grid_ready||sn_grid_T!=T||sn_grid_B!=B||sn_grid_width!=sn_width) {
        sn_grid_ready=0;
        for(int i=0;i<=N_KAPPA_F;i++) {
            sn_x[i]=log(NUMIN)+i*log(NUMAX/NUMIN)/N_KAPPA_F;
            double f=kappa_synch_num_F_direct(T,B,exp(sn_x[i]));
            if(f<=0.||!isfinite(f))synch_num_fail("unrepresentable F grid: use direct path",f);
            sn_y[i]=log(f);
        }
        double maxerr=0.;
        for(int i=0;i<N_KAPPA_F;i++) {
            double f=kappa_synch_num_F_direct(T,B,exp(.5*(sn_x[i]+sn_x[i+1])));
            double err=fabs(exp(.5*(sn_y[i]+sn_y[i+1]))/f-1.);
            if(!isfinite(err)||err>SYNCH_F_INTERP_TOL)synch_num_fail("F interpolation accuracy",err);
            maxerr=fmax(maxerr,err);
        }
        sn_grid_T=T;sn_grid_B=B;sn_grid_width=sn_width;sn_grid_ready=1;
        fprintf(stderr,"%s T=%.17g width=%.17g B=%.17g gamma_cut=%.17g upper=%.17g photon_cut=%.17g electron_order=%d angle_order=%d F_intervals=%d midpoint_max=%.17g\n",
                SYNCH_NUM_ID,T,sn_width,B,(double)KAPPA_GAMMA_CUT,100.*KAPPA_GAMMA_CUT,(double)nu_cutoff,KAPPA_SYNCH_ELECTRON_ORDER,KAPPA_SYNCH_ANGLE_ORDER,N_KAPPA_F,maxerr);
    }
    for(int i=0;i<=N_KAPPA_F;i++){x[i]=sn_x[i];y[i]=sn_y[i];}
    return N_KAPPA_F;
}
static double synch_num_F(double T,double B,double nu) {
    double f;
    if(!isfinite(nu)||nu<=0.)synch_num_fail("invalid F frequency",nu);
    if(sn_grid_ready&&T==sn_grid_T&&B==sn_grid_B&&nu>=NUMIN&&nu<=NUMAX) {
        double v=(log(nu)-log(NUMIN))*N_KAPPA_F/log(NUMAX/NUMIN);
        int i=(int)floor(v);if(i>=N_KAPPA_F)i=N_KAPPA_F-1;
        double a=(log(nu)-sn_x[i])/(sn_x[i+1]-sn_x[i]);
        f=exp((1.-a)*sn_y[i]+a*sn_y[i+1]);
    } else f=kappa_synch_num_F_direct(T,B,nu);
    return f*exp(-nu/nu_cutoff);
}
