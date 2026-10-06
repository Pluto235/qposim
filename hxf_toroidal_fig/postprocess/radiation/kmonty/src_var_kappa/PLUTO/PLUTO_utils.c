
#include "decs.h"
#include "PLUTO_model.h"
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* Harm globals */
extern double ****econ;
extern double ****ecov;
extern double ***bcon;
extern double ***bcov;
extern double ***ucon;
extern double ***ucov;
extern double ****p;
extern struct of_geom **geom;
extern double **ne;
extern double **thetae;
extern double **b;
double weight_min;

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
void coord(int i, int j,int k ,double *X);
void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
                    double *sigma, double *beta, double Ucon[NDIM],
                    double Bcon[NDIM], int *ACCZONE, double *dx_local,
                    int *igrid_c);

/** HARM utilities **/

/********************************************************************

   Interpolation routines

********************************************************************/

double interp_scalar(double ***var, int i, int j, int k, double coeff[4]) {

    double interp;
    int ip1, jp1, kp1;
    double b1, b2, b3, del[NDIM];

    del[1] = coeff[1];
    del[2] = coeff[2];
    del[3] = coeff[3];
    if (del[1] > 1+1e-4 || del[2] > 1+1e-4 || del[3] > 1+1e-4 || del[1] < -1e-4 || del[2] < -1e-4 || del[3] < -1e-4){
        printf("Over boundary: i=%d j=%d k=%d del[1]=%.4e del[2]=%.4e del[3]=%.4e\n",i,j,k,del[1],del[2],del[3]);
        return -1;
        }

    ip1 = i + 1;
    jp1 = j + 1;
    kp1 = k + 1;

    b1 = 1. - del[1];
    b2 = 1. - del[2];
    b3 = 1. - del[3];
    //turn off interpolation outside the boundary
    if (k==-1)return var[0][0][0];
    else if (k==N3)return var[0][0][N3-1];

    interp = var[i][j][k] * b1 * b2 + var[i][jp1][k] * b1 * del[2] +
             var[ip1][j][k] * del[1] * b2 + var[ip1][jp1][k] * del[1] * del[2];

    /* Now interpolate above in x3 */
    interp = b3 * interp + del[3] * (var[i][j][kp1] * b1 * b2 +
                                     var[i][jp1][kp1] * b1 * del[2] +
                                     var[ip1][j][kp1] * del[1] * b2 +
                                     var[ip1][jp1][kp1] * del[1] * del[2]);

    return interp;
}

double lnu_min, lnu_max, dlnu;

static void init_linear_interp_weight() {

    lnu_min = log(NUMIN);
    lnu_max = log(NUMAX);
    dlnu = (lnu_max - lnu_min) / (N_ESAMP);
}

static double linear_interp_weight(double nu) {
    int i;
    double di, lnu;

    lnu = log(nu);

    di = (lnu - lnu_min) / dlnu;
    i = (int)(di + 2.5) - 2;
    di = di - i;

    return exp((1. - di) * wgt[i] + di * wgt[i + 1]); //* exp(-nu/nu_cutoff);
}

/***********************************************************************************

   End interpolation routines

***********************************************************************************/

#define JCST (M_SQRT2 * EE * EE * EE / (27 * ME * CL * CL))
void init_weight_table(void) {

    int i, j, k, l, lstart, lend, myid, nthreads;
    double Ne, Thetae, B, K2;
    double sum[N_ESAMP + 1], nu[N_ESAMP + 1];
    double fac_th, fac_nth, sfac;
    double Ucon[NDIM], Bcon[NDIM];

    /*      Set up interpolation */
    init_linear_interp_weight();

#pragma omp parallel for schedule(static) private(i)
    for (i = 0; i <= N_ESAMP; i++) {
        sum[i] = 0.;
        sum[i] = 0.;
        nu[i] = exp(i * dlnu + lnu_min);
    }

    
#pragma omp parallel private(i, j, k, Thetae, K2, Ne, B, fac_th, fac_nth, l,   \
                             lstart, lend, myid, nthreads, Ucon, Bcon)
    {
        nthreads = omp_get_num_threads();
        myid = omp_get_thread_num();
        lstart = myid * (N_ESAMP / nthreads);
        lend = (myid + 1) * (N_ESAMP / nthreads);
        if (myid == nthreads - 1)
            lend = N_ESAMP + 1;
        int ACCZONE;
        for (i = 0; i < N1; i++)
            for (j = 0; j < N2; j++) {
                for (k = 0; k < N3; k++) {
					sfac = dx[1][i]* dx[2][j] * dx[3][k] * L_unit * L_unit * L_unit;
                    get_fluid_zone(i, j, k, &Ne, &Thetae, &B, Ucon, Bcon,&ACCZONE);
                    if (Ne == 0. || Thetae < THETAE_MIN)
                        continue;

                        // thermal
#if (THERMAL)
                    K2 = K2_eval(Thetae);
                    fac_th = (JCST * Ne * B * Thetae * Thetae / K2) * sfac *
                             geom[i][j].g;
#elif (KAPPA || POWERLAW)
                    fac_nth = (EE * EE * EE / (2. * M_PI * ME * CL * CL)) * Ne *
                              B * sfac * geom[i][j].g;
#elif (MIXED)
                    K2 = K2_eval(Thetae);
                    fac_th = (JCST * Ne * B * Thetae * Thetae / K2) * sfac *
                             geom[i][j].g;
                    fac_nth = (EE * EE * EE / (2. * M_PI * ME * CL * CL)) * Ne *
                              B * sfac * geom[i][j].g;
#endif
                    for (l = lstart; l < lend; l++) {
#if (THERMAL)
                        sum[l] += fac_th * F_eval_th(Thetae, B, nu[l]);
#elif (KAPPA || POWERLAW)
                        sum[l] += fac_nth * F_eval_nth(Thetae, B, nu[l]);
#elif (MIXED)

                        if (ACCZONE) {
                            sum[l] += perct_thermal * fac_th *
                                          F_eval_th(Thetae, B, nu[l]) +
                                      (1 - perct_thermal) * fac_nth *
                                          F_eval_nth(Thetae, B, nu[l]);
                        } else
                            sum[l] += fac_th * F_eval_th(Thetae, B, nu[l]);
#endif
						if (isnan(sum[l]))
                            fprintf(stderr,
                                    "fac_th %e K2 %e Thetae %e B %e Ne %e "
                                    "nu[l] %e F_eval %e\n",
                                    fac_th, K2, Thetae, B, Ne, nu[l],
                                    F_eval_th(Thetae, B, nu[l]));
                    }
                }
            }
#pragma omp barrier
    }
    double weight_max=0.;
#pragma omp parallel for schedule(static) private(i)
    for (i = 0; i <= N_ESAMP; i++) {
        wgt[i] = log(sum[i]/(HPL * Ns));
        if (wgt[i] > weight_max) weight_max = wgt[i];
    }
    weight_min=weight_max/WEIGHT_MIN;
    /* Write weights to file for debugging */
    // FILE *fp = fopen("weight.txt", "w");
    // if (fp != NULL) {
    //     for (i = 0; i <= N_ESAMP; i++) {
    //         fprintf(fp, "%e %e\n", nu[i], wgt[i]);
    //     }
    //     fclose(fp);
    //     fprintf(stderr, "Weight table written to weight.txt\n");
    // }
    // exit(1);
    fflush(stderr);
   
    return;
}

#undef JCST

double lb_min, dlb;
double nint[NINT + 1];
double dndlnu_max[NINT + 1];

void init_nint_table(void) {
    int i, j;
    double Bmag, dn = 0;
    static int firstc = 1;

    if (firstc) {
        lb_min = log(BTHSQMIN);
        dlb = log(BTHSQMAX / BTHSQMIN) / NINT;
        firstc = 0;
    }

    for (i = 0; i <= NINT; i++) {
        nint[i] = 0.;
        Bmag = exp(i * dlb + lb_min);
        dndlnu_max[i] = 0.;
        for (j = 0; j < N_ESAMP; j++) {
#if (THERMAL)
            dn = F_eval_th(1., Bmag, exp(j * dlnu + lnu_min)) /
                 (exp(wgt[j]) + 1.e-100);
#elif (KAPPA || POWERLAW)
            dn = F_eval_nth(1., Bmag, exp(j * dlnu + lnu_min)) /
                 (exp(wgt[j]) + 1.e-100);

#endif
            if (dn > dndlnu_max[i]) dndlnu_max[i] = dn;
            nint[i] += dlnu * dn;
        }
#if (THERMAL)
        double fac_th = M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL);

        nint[i] *=L_unit * L_unit * L_unit * fac_th * 1. / HPL;
#elif (KAPPA || POWERLAW)
        double fac_nth = (EE * EE * EE / (2. * M_PI * ME * CL * CL));
        nint[i] *=L_unit * L_unit * L_unit * fac_nth * 1. / HPL;
#endif
        nint[i] = log(nint[i]);
        dndlnu_max[i] = log(dndlnu_max[i]);
    }

    return;
}
#undef JCST

static void init_zone(int i, int j, int k, double *nz, double *dnmax) {

    int l;
    double Ne, Thetae, Bmag, lbth, fac_th, fac_nth;
    double dl, dn, ninterp, K2;
    double Ucon[NDIM], Bcon[NDIM];
    int ACCZONE = 0;
    double r,X[NDIM];
    coord(i,j,k,X);
    r=sqrt(X[1]*X[1]+X[2]*X[2]);
    if (r>RMAX){//stop radius
        *nz = 0.;
        *dnmax = 0.;
        return;
    }

    get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon, &ACCZONE);

    if (Ne == 0. || Thetae < THETAE_MIN || Bmag < 1e-20) {
        *nz = 0.;
        *dnmax = 0.;
        return;
    }
#if (POWERLAW)
    lbth = log(Bmag);
#else
    lbth = log(Bmag * Thetae * Thetae);
#endif
    lb_min = log(BTHSQMIN);
    dlb = log(BTHSQMAX / BTHSQMIN) / NINT;

    dl = (lbth - lb_min) / dlb;
    l = (int)dl;
    dl = dl - l;
    if (l < 0) {
        *dnmax = 0.;
        *nz = 0.;
		//fprintf(stderr,
        //        "l smaller than zero l %d lbth %e dlb %e Bmag %e Thetae %e\n",
        //        l, lbth, dlb, Bmag, Thetae);
        return;
    } 
    else if (l >= NINT || MIXED || POWERLAW) {
		ninterp = 0.;
        *dnmax = 0.;
        for (l = 0; l <= N_ESAMP; l++) { 
#if (KAPPA || POWERLAW)
            dn = F_eval_nth(Thetae, Bmag, exp(l * dlnu + lnu_min)) / exp(wgt[l]);
#elif (THERMAL)
            dn = F_eval_th(Thetae, Bmag, exp(l * dlnu + lnu_min)) / exp(wgt[l]);
#elif (MIXED)
            K2 = K2_eval(Thetae);
            if (K2 == 0.) {
                *nz = 0.;
                *dnmax = 0.;
                return;
            }

            if (ACCZONE) {
                dn = perct_thermal *
                         F_eval_th(Thetae, Bmag, exp(l * dlnu + lnu_min)) /
                         (exp(wgt[l] + 1e-100)) +
                     (27. / (M_SQRT2 * 2. * M_PI)) * (K2 / (Thetae * Thetae)) *
                         (1 - perct_thermal) *
                         F_eval_nth(Thetae, Bmag, exp(l * dlnu + lnu_min)) /
                         (exp(wgt[l] + 1e-100));
            } else
                dn = F_eval_th(Thetae, Bmag, exp(l * dlnu + lnu_min)) /
                     (exp(wgt[l]));
#endif

            if (dn > *dnmax) {
                *dnmax = dn;
            }
            ninterp += dlnu * dn;
        }
#if (THERMAL)
        fac_th = M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL);
        ninterp *= dx[1][i] * dx[2][j] * dx[3][k] * L_unit * L_unit * L_unit *
                   fac_th * 1. / HPL;
        //ninterp *= L_unit * L_unit * L_unit *fac_th * 1. / HPL;
#elif (KAPPA || POWERLAW)
        fac_nth = (EE * EE * EE / (2. * M_PI * ME * CL * CL));
        ninterp *= dx[1][i] * dx[2][j] * dx[3][k] * L_unit * L_unit * L_unit *
                   fac_nth * 1. / HPL;

#elif (MIXED)
        double fac_mixed = M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL);

        ninterp *= dx[1][i] * dx[2][j] * dx[3][k] * L_unit * L_unit * L_unit *
                   fac_mixed * 1. / HPL;

#endif
    } 
    
    else {
        if (isinf(nint[l]) || isinf(nint[l + 1])) {
            ninterp = 0.;
            *dnmax = 0.;
        } 
        else {
            ninterp=(exp((1. - dl) * nint[l] + dl * nint[l + 1]))*dx[1][i]*dx[2][j]*dx[3][k];
            *dnmax = exp((1. - dl) * dndlnu_max[l] + dl * dndlnu_max[l + 1]);
        }
    }

#if (THERMAL)
    K2 = K2_eval(Thetae);
    if (K2 == 0.) {
        *nz = 0.;
        *dnmax = 0.;
        return;
    }

    fac_th = Ne * Bmag * Thetae * Thetae / K2;
    *nz = geom[i][j].g * fac_th * ninterp;
#elif (KAPPA || POWERLAW)
    fac_nth = Ne * Bmag;
    *nz = geom[i][j].g * fac_nth * ninterp;

#elif (MIXED)
    double fac_mixed = Ne * Bmag * Thetae * Thetae / K2;
    *nz = geom[i][j].g * fac_mixed * ninterp;
#endif

    if (*nz > Ns * log(NUMAX / NUMIN)) {
        // fprintf(stderr,"\n %e %e %e %e %e %e %e %e\n", p[KRHO][i][j][k],
        //         p[UU][i][j][k], p[U1][i][j][k], p[U2][i][j][k],
        //         p[U3][i][j][k], p[B1][i][j][k], p[B2][i][j][k],
        //         p[B3][i][j][k]);
        fprintf(stderr,
                "Something very wrong in zone %d %d %d: \n Ne=%g B=%g  Thetae=%g  K2=%g  "
                "ninterp=%g *nz=%g\n\n",
                i, j,k,Ne, Bmag ,Thetae, K2, ninterp,*nz);
        *nz = 0.;
        *dnmax = 0.;
    }
    return;
}

int zone_flag;
int get_zone(int *i, int *j, int *k, double *dnmax) {
    /* Return the next zone and the number of superphotons that need to be
     * generated in it.  */
    int in2gen;
    double n2gen;
    static int zi = 0;
    static int zj = 0;
    static int zk = -1;
    zone_flag = 1;
    zk++;
    if (zk >= N3) {
        zk = 0;
        zj++;
        if (zj >= N2) {
            zj = 0;
            zi++;
            if (zi >= N1) {
                *i = N1;
                return 1;
            }
        }
    }
    if (zi >= N1) {
        *i = N1;
        return 1;
    }
    init_zone(zi, zj, zk, &n2gen, dnmax);
    if (fmod(n2gen, 1.) > monty_rand()) {
        in2gen = (int)n2gen + 1;
    } else {
        in2gen = (int)n2gen;
    }

    *i = zi;
    *j = zj;
    *k = zk;

    // we only need a fraction of these for every processor

    return in2gen;
}
void sample_zone_photon(int i, int j, int k, double dnmax,
                        struct of_photon *ph) {
    /* Set all initial superphoton attributes */

    int l;
    double K_tetrad[NDIM], tmpK[NDIM], E, Nln;
    double nu, th, cth, sth, phi, sphi, cphi, jmax, weight;
    double Ne, Thetae, Bmag, Ucon[NDIM], Bcon[NDIM], bhat[NDIM];
    static double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
    double sigma, beta, dx_local;
    int igrid = i;

    coord(i, j,k, ph->X);

    Nln = lnu_max - lnu_min;
    int ACCZONE = 0;
    get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, &sigma, &beta, Ucon, Bcon, &ACCZONE, &dx_local, &igrid);
	
    /* Sample from superphoton distribution in current simulation zone */
	
    do {
        nu = exp(monty_rand() * Nln + lnu_min);
        weight = linear_interp_weight(nu);
    } while (monty_rand() >
             (F_eval(Thetae, Bmag, nu, ACCZONE) / weight) / dnmax);

    ph->w = weight;
    jmax = jnu_synch(nu, Ne, Thetae, Bmag, M_PI / 2., ACCZONE, sigma, beta);
	
    do {
        cth = 2. * monty_rand() - 1.;
        th = acos(cth);

    } while (monty_rand() >
             jnu_synch(nu, Ne, Thetae, Bmag, th, ACCZONE, sigma, beta) / jmax);
	
    sth = sqrt(1. - cth * cth);
    phi = 2. * M_PI * monty_rand();
    cphi = cos(phi);
    sphi = sin(phi);

    E = nu * HPL / (ME * CL * CL);
    K_tetrad[0] = E;
    K_tetrad[1] = E * cth;
    K_tetrad[2] = E * cphi * sth;
    K_tetrad[3] = E * sphi * sth;

 
    if (zone_flag) { /* first photon created in this zone, so make the tetrad */
        if (Bmag > 0.) {
            for (l = 0; l < NDIM; l++)
                bhat[l] = Bcon[l] * B_unit / Bmag;
        } else {
            for (l = 1; l < NDIM; l++)
                bhat[l] = 0.;
            bhat[1] = 1.;
        }
        make_tetrad(ph->X, Ucon, bhat, geom[i][j].gcov, Econ, Ecov);
        zone_flag = 0;
    }

    tetrad_to_coordinate(Econ, K_tetrad, ph->K);

    K_tetrad[0] *= -1.;
    tetrad_to_coordinate(Ecov, K_tetrad, tmpK);

    ph->E = ph->E0 = ph->E0s = -tmpK[0];
    ph->L = tmpK[3];
    ph->tau_scatt = 0.;
    ph->tau_abs = 0.;
    ph->X1i = ph->X[1];
    ph->X2i = ph->X[2];
    ph->X3i = ph->X[3];
    ph->nscatt = 0;
    ph->ne0 = Ne;
    ph->b0 = Bmag;
    ph->thetae0 = Thetae;
	
    return;
}

int binary_search(double arr[],double arr2[], int low,int high, double target) {
	if (target<arr[0])return -1;
	if (target>arr2[high])return high;
    while (low <= high) {
        int mid = low+(high - low)/2;
        if ((arr[mid]+arr2[mid])/2 < target) low = mid + 1;
        else if ((arr[mid]+arr2[mid])/2 > target) high = mid - 1;
        else return (mid-1);
    }
    if(target>(arr[high]+arr2[high])/2) return high;
    else return (high-1);//smaller than middle, move to the left grid
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]) {
    int low,high;
	
  	if (fabs(dx1max-dx1min)<1.e-4) *i=(int)floor((X[1]-x1l[0])/dx1min-0.5);
    else {// predict coordinate
      	low=(int)floor((X[1]-x1l[0])/dx1max-0.5);
      	high=(int)floor((X[1]-x1l[0])/dx1min+0.5);
      	*i = binary_search(x1l,x1r,low,MIN(high,N1-1),X[1]);
        }
  	if (fabs(dx2max-dx2min)<1.e-4) *j=(int)floor((X[2]-x2l[0])/dx2min-0.5);
    else {
      	low=(int)floor((X[2]-x2l[0])/dx2max-0.5);
      	high=(int)floor((X[2]-x2l[0])/dx2min+0.5);
      	*j = binary_search(x2l,x2r,low,MIN(high,N2-1),X[2]);
        }
  	if (fabs(dx3max-dx3min)<1.e-4) *k=(int)floor((X[3]-x3l[0])/dx3min-0.5);
    else {
      	low=(int)floor((X[3]-x3l[0])/dx3max-0.5);
      	high=(int)floor((X[3]-x3l[0])/dx3min+0.5);
      	*k = binary_search(x3l,x3r,low,MIN(high,N3-1),X[3]);
        }
	
    if (*i < 0) {
        *i = 0;
        del[1] = 0.;
    } else if (*i > N1 - 2) {
        *i = N1 - 2;
        del[1] = 1.;
    } else {
        del[1] = (X[1]-(x1l[*i]+x1r[*i])/2)/((x1r[*i]-x1l[*i]+x1r[*i+1]-x1l[*i+1])/2);
    }

    if (*j < 0) {
        *j = 0;
        del[2] = 0.;
    } else if (*j > N2 - 2) {
        *j = N2 - 2;
        del[2] = 1.;
    } else {
        del[2] = (X[2]-(x2l[*j]+x2r[*j])/2) /((x2r[*j]-x2l[*j]+x2r[*j+1]-x2l[*j+1])/2);
    }

    if (*k < 0) {
        *k = -1;//xf: as a flag
        del[3] = 0.;
    } else if (*k > N3 - 2) {
        *k = N3 ;//xf: as a flag
        del[3] = 1.;
    } else {
        del[3] = (X[3]-(x3l[*k]+x3r[*k])/2) /((x3r[*k]-x3l[*k]+x3r[*k+1]-x3l[*k+1])/2);
    }
}


void coord(int i, int j,int k, double *X) {
    /* returns zone-centered values for coordinates */
    X[0] = startx[0];
    X[1] = (x1l[i]+x1r[i])/2.;
    X[2] = (x2l[j]+x2r[j])/2.;
    X[3] = (x3l[k]+x3r[k])/2.;
    return;
}

void set_units(char *munitstr) {
    /** from this, calculate units of length, time, mass,
        and derivative units **/
    L_unit = GNEWT * MBH / (CL * CL);
    T_unit = L_unit / CL;

    RHO_unit = M_unit / pow(L_unit, 3);
    U_unit = RHO_unit * CL * CL;
    B_unit = CL * sqrt(4. * M_PI * RHO_unit);
    Ne_unit = RHO_unit / (MP + ME);

    max_tau_scatt = 1e-5 * (6. * L_unit) * RHO_unit * 0.4;
    // Get the rank of the process
    int world_rank;
#if MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) 
#endif
	{
        fprintf(stderr, "\nUNITS\n");
		fprintf(stderr, "MBH \t\t= %g Msun\n", MBH/MSUN);
		fprintf(stderr, "DISTANCE \t= %g kpc\n", SOURCE_DIST/KPC2CM);
        fprintf(stderr, "L_unit,T_unit,M_unit: %g cm %g s %g g\n", L_unit, T_unit,M_unit);
        fprintf(stderr, "RHO_unit,U_unit,B_unit: %g g/cm^3 %g g/(cm s^2) %g Gauss\n",RHO_unit,U_unit,B_unit);
        fprintf(stderr, "Ne: %g \n", Ne_unit);
        fprintf(stderr, "max_tau_scatt: %g\n", max_tau_scatt);
    }
}

/* set up all grid functions */
void init_geometry() {
    int i, j;
    double X[NDIM]={0,0,0,0};

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {

            /* zone-centered */
            //coord(i, j, X);

            gcov_func(X, geom[i][j].gcov);

            geom[i][j].g = gdet_func(geom[i][j].gcov, X);

            gcon_func(geom[i][j].gcov, geom[i][j].gcon);
        }
    }

    /* done! */
}

/*

   return solid angle between points x2i, x2f
   and over all x3.

 */

double dOmega_func(double x2i, double x2f, double x3i, double x3f) {
    double dO;

    dO = (x3f - x3i) * (-cos(M_PI * x2f) + cos(M_PI * x2i));

    return (dO);
}

static void *malloc_rank1(int n1, int size) {
    void *A;

    if ((A = (void *)malloc(n1 * size)) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank1\n");
        exit(123);
    }

    return A;
}

static void **malloc_rank2(int n1, int n2, int size) {

    void **A;
    int i;

    if ((A = (void **)malloc(n1 * sizeof(void *))) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank2\n");
        exit(124);
    }

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, size);
    }

    return A;
}

static double **malloc_rank2_cont(int n1, int n2) {

    double **A;
    double *space;
    int i;

    space = malloc_rank1(n1 * n2, sizeof(double));

    A = malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++)
        A[i] = &(space[i * n2]);

    return A;
}

double ***malloc_rank3(int n1, int n2, int n3) {
    double ***A;
    double *space;
    int i, j;

    space = malloc_rank1(n1 * n2 * n3, sizeof(double));

    A = malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, sizeof(double *));
        for (j = 0; j < n2; j++) {
            A[i][j] = &(space[n3 * (j + n2 * i)]);
        }
    }

    return A;
}

void init_storage(void) {
    int i;

    p = (double ****)malloc_rank1(NPRIM, sizeof(double *));
    for (i = 0; i < NPRIM; i++)
        p[i] = malloc_rank3(N1, N2, N3);

    geom = (struct of_geom **)malloc_rank2(N1, N2, sizeof(struct of_geom));

    int j, k;
    shared_spect = (struct of_spectrum ***)malloc(
        N_THBINS * sizeof(struct of_spectrum **));
    for (i = 0; i < N_THBINS; i++) {
        shared_spect[i] = (struct of_spectrum **)malloc(
            N_PHIBINS * sizeof(struct of_spectrum *));
        for (j = 0; j < N_PHIBINS; j++) {
            shared_spect[i][j] = (struct of_spectrum *)malloc(
                N_EBINS * sizeof(struct of_spectrum));
        }
    }

    return;
}
