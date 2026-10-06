/*

   sphere model specification routines
   //model from B. Ryan and C. Gammie

 */

#include "decs.h"
#define global
#include "PLUTO_model.h"
#undef global
#define MAX(a, b) ((a) >(b) ? (a) : (b))


struct of_spectrum ***spect;
extern struct of_spectrum ***shared_spect;

double *dx[NDIM];
int N1, N2, N3;
double gam;

#pragma omp threadprivate(spect)

/*

   encapsulates all initialization routines

 */

void init_model(char *args[]) {
    /* find dimensional quantities from black hole
       mass and its accretion rate */
    set_units();

    int world_rank;
#if MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
#endif
        fprintf(stderr, "READING DATA\n");

    init_pluto_data(args[2]); /* read in pluto simulation data */

#if MPI
    if (world_rank == 0)
#endif
		fprintf(stderr, "bias_norm: %g\n", bias_norm);

#pragma omp parallel
    {
        int i, j, k;
        spect = (struct of_spectrum ***)malloc(N_THBINS *
                                               sizeof(struct of_spectrum **));
        for (i = 0; i < N_THBINS; i++) {
            spect[i] = (struct of_spectrum **)malloc(
                N_PHIBINS * sizeof(struct of_spectrum *));
            for (j = 0; j < N_PHIBINS; j++) {
                spect[i][j] = (struct of_spectrum *)malloc(
                    N_EBINS * sizeof(struct of_spectrum));
            }
        }
    }

    init_geometry();
    if (world_rank == 0)
        fprintf(stderr, "init geom done\n");
    /* make look-up table for hot cross sections */
    init_hotcross();
    if (world_rank == 0)
        fprintf(stderr, "hot cross done\n");
    /* make table for solid angle integrated emissivity and K2 */
    init_emiss_tables();
    if (world_rank == 0)
        fprintf(stderr, "emiss tables done\n");

    /* make table for superphoton weights */
    init_weight_table();
    if (world_rank == 0)
        fprintf(stderr, "weight tables done\n");

    /* make table for quick evaluation of ns_zone */
    init_nint_table();
	if (world_rank == 0)
        fprintf(stderr, "ns tables done\n");
}

/*
   make super photon
 */

int n2gen = -1;
double dnmax;
int zone_i, zone_j, zone_k;

void make_super_photon(struct of_photon *ph, int *quit_flag) {
    // fprintf(stderr,"making some superphotons...\n");

    while (n2gen <= 0) {
        // fprintf(stderr,"n2gen is %d %d\n",zone_i,n2gen);
        /*dnmax=number of superphotons needed to be crated in a given zone*/
        n2gen = get_zone(&zone_i, &zone_j, &zone_k, &dnmax);

        //fprintf(stdout,"i=%d j=%d k=%d n2gen=%d dnmax=%g\n",
		//		zone_i,zone_j,zone_k,n2gen,dnmax);
    }

    n2gen--;

    if (zone_i == N1)
        *quit_flag = 1;
    else
        *quit_flag = 0;

    if (*quit_flag != 1) {
        /* Initialize the superphoton energy, direction, weight, etc. */
		//fprintf(stdout,"i=%d j=%d k=%d\n",zone_i,zone_j,zone_k);
        sample_zone_photon(zone_i, zone_j, zone_k, dnmax, ph);
        // fprintf(stderr,"X= %g %g %g K=%g %g %g\n",
        //         ph->X[1],ph->X[2],ph->X[3],ph->K[1],ph->K[2],ph->K[3]);
	#if MPI
		int world_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
		if (world_rank == 0)
	#endif
		{
		if (zone_i%40==0 && zone_j%100==0)
			fprintf(stderr,"[progress] %.2lf \%\n",
				    (float)(zone_i*N2+zone_j)/(float)(N1*N2)*100.);
    	}
    }
    return;
}

/*

   produces a bias (> 1) for probability of Compton scattering
   as a function of local temperature

 */

double bias_func(double Te, double w) {
    double bias;
    bias = 100. * Te * Te / bias_norm;
#pragma omp atomic
    bias /= max_tau_scatt;

    if (bias < 1.) bias = 1.;
	return bias;
	//else if (bias<100) return bias;
    //else return sqrt(100.*bias);
}

void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
                    double *sigma, double *beta, double Ucon[NDIM],
                    double Bcon[NDIM], int *ACCZONE, double *dx_local,
                    int *igrid_c) {
    double Ucov[NDIM], Bcov[NDIM];
    double Bp[NDIM], Vcon[NDIM], Vfac, V[NDIM],lor,VdotV, UdotBp;
    double X[NDIM];
    coord(i, j,k, X);
    // int ti=i,tj=j,tk=k;
    *Ne = p[KRHO][i][j][k] * Ne_unit;

    Bp[1] = p[B1][i][j][k];
    Bp[2] = p[B2][i][j][k];
    Bp[3] = p[B3][i][j][k];

    V[1] = p[U1][i][j][k];
    V[2] = p[U2][i][j][k];
    V[3] = p[U3][i][j][k];

    lor=1./sqrt(1.-V[1]*V[1]-V[2]*V[2]-V[3]*V[3]);
    Vcon[1]=lor*V[1];
    Vcon[2]=lor*V[2];
    Vcon[3]=lor*V[3];

    // Get Ucov
    VdotV = 0.;
    for (int l = 1; l < NDIM; l++)
        for (int m = 1; m < NDIM; m++)
            VdotV += geom[i][j].gcov[l][m] * Vcon[l] * Vcon[m];
    Vfac = sqrt(-1. / geom[i][j].gcon[0][0] * (1. + fabs(VdotV)));
    Ucon[0] = -Vfac * geom[i][j].gcon[0][0];
    for (int l = 1; l < NDIM; l++)
        Ucon[l] = Vcon[l] - Vfac * geom[i][j].gcon[0][l];
    lower(Ucon, geom[i][j].gcov, Ucov);
    // Get B and Bcov
    UdotBp = 0.;
    for (int l = 1; l < NDIM; l++)
        UdotBp += Ucov[l] * Bp[l];
    Bcon[0] = UdotBp;
    for (int l = 1; l < NDIM; l++)
        Bcon[l] = (Bp[l] + Ucon[l] * UdotBp) / Ucon[0];

    lower(Bcon, geom[i][j].gcov, Bcov);

    *B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
              Bcon[3] * Bcov[3]) * B_unit +1e-30;

	double Bsq = pow(*B/B_unit,2.);
	double b2 = pow(p[UU][i][j][k] * (gam - 1.) / (0.5 * (Bsq + SMALL)), 2.);
	double trat = R_HIGH * b2 / (1. + b2) + R_LOW / (1. + b2);
	Thetae_unit = (gam - 1.) * ((MP / ME)) / (trat + 1);

	*Thetae = (p[UU][i][j][k] / p[KRHO][i][j][k]) * Thetae_unit;

    // Compute sigma and beta parameters for kappa distribution
    double rho = p[KRHO][i][j][k];
    double uu = p[UU][i][j][k];
    *sigma = Bsq / rho;
    *beta = uu * (gam - 1.) / (0.5 * (Bsq + SMALL));

    // i=ti,j=tj,k=tk;
    // if(i==245 &&j==218 && k==409){
    //     fprintf(stderr,"B=%g b2=%g Thetae=%g thetae_unit=%g trat=%g\n",
    //            *B,b2,*Thetae,Thetae_unit,trat);
    //     fprintf(stderr,"\n %e %e %e %e %e %e %e %e\n", p[KRHO][i][j][k],
    //            p[UU][i][j][k], p[U1][i][j][k], p[U2][i][j][k],
    //            p[U3][i][j][k], p[B1][i][j][k], p[B2][i][j][k],
    //            p[B3][i][j][k]);
    //     fprintf(stderr,"U[0]=%g U[1]=%g U[2]=%g U[3]=%g\n",
    //            Ucon[0],Ucon[1],Ucon[2],Ucon[3]);
    //     fprintf(stderr,"Bcon[0]=%g Bcon[1]=%g Bcon[2]=%g Bcon[3]=%g\n",
    //            Bcon[0],Bcon[1],Bcon[2],Bcon[3]);
    //     //exit(1);
    // }

    //sig = pow(*B / B_unit, 2) / (*Ne / Ne_unit);
    //if (sig > 1.) {
    //    *Thetae = SMALL; // *Ne = 1.e-10*Ne_unit;
    //    *Ne = 0;
        // fprintf(stderr, "yup its high!\n");
    //}
}

int get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
                     double *Thetae, double *B, double *sigma, double *beta,
                     double Ucon[NDIM], double Ucov[NDIM], double Bcon[NDIM],
                     double Bcov[NDIM], int *ACCZONE, double *dx_local,
                     int *igrid_c) {
    int i, j, k;
    double rho, uu;
    double del[NDIM];
    double Bp[NDIM], Vcon[NDIM],V[NDIM],lor, Vfac, VdotV, UdotBp;
    double gcon[NDIM][NDIM], coeff[4];
    double interp_scalar(double ***var, int i, int j,int k, double del[4]);

    if (X[1] < startx[1] || X[1] > stopx[1] || X[2] < startx[2] ||
        X[2] > stopx[2] || X[3] < startx[3] || X[3] > stopx[3]) {
        *Ne = 0;
        return 0;
    }

    Xtoijk(X, &i, &j, &k, del);
    // int ti=i,tj=j,tk=k;

    rho = interp_scalar(p[KRHO], i, j, k, del);
    uu = interp_scalar(p[UU], i, j,k, del);

	*Ne = rho * Ne_unit + SMALL;

    Bp[1] = interp_scalar(p[B1], i, j,k, del);
    Bp[2] = interp_scalar(p[B2], i, j,k, del);
    Bp[3] = interp_scalar(p[B3], i, j,k, del);

    V[1] = interp_scalar(p[U1], i, j,k, del);
    V[2] = interp_scalar(p[U2], i, j,k, del);
    V[3] = interp_scalar(p[U3], i, j,k, del);

    lor=1./sqrt(1.-V[1]*V[1]-V[2]*V[2]-V[3]*V[3]);
    Vcon[1]=lor*V[1];
    Vcon[2]=lor*V[2];
    Vcon[3]=lor*V[3];

    gcov_func(X, gcov);
    gcon_func(gcov, gcon);

    // Get Ucov
    VdotV = 0.;
    for (int i = 1; i < NDIM; i++)
        for (int j = 1; j < NDIM; j++)
            VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
    Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
    Ucon[0] = -Vfac * gcon[0][0];
    for (int i = 1; i < NDIM; i++)
        Ucon[i] = Vcon[i] - Vfac * gcon[0][i];
    lower(Ucon, gcov, Ucov);

    // Get B and Bcov
    UdotBp = 0.;
    for (int i = 1; i < NDIM; i++)
        UdotBp += Ucov[i] * Bp[i];
    Bcon[0] = UdotBp;
    for (int i = 1; i < NDIM; i++)
        Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];

    lower(Bcon, gcov, Bcov);

    *B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
              Bcon[3] * Bcov[3]) *B_unit +1e-30;

	double Bsq = fabs(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
                      Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]);

    *beta = uu * (gam - 1.) / (0.5 *Bsq + SMALL);

    double b2 = pow(*beta, 2.);

    *sigma = Bsq / (rho);

    double trat = R_HIGH * b2 / (1. + b2) + R_LOW / (1. + b2);
	//double trat=3.;

    Thetae_unit = (gam - 1.) * ((MP / ME)) / (trat + 1);

    *Thetae = (uu / rho) * Thetae_unit;

//    i=ti,j=tj,k=tk;
//    if(i==201 &&j==250 && k==598){
//        fprintf(stderr,"beta=%g sigma=%g thetae=%g thetae_unit=%g\n",
//               *beta,*sigma,*Thetae,Thetae_unit);
//        fprintf(stderr,"\n %e %e %e %e %e %e %e %e\n", p[KRHO][i][j][k],
//               p[UU][i][j][k], p[U1][i][j][k], p[U2][i][j][k],
//               p[U3][i][j][k], p[B1][i][j][k], p[B2][i][j][k],
//               p[B3][i][j][k]);
//        //exit(1);
//    }

    //if (sig > 1.) {
    //    *Thetae = SMALL; // *Ne = 1.e-10*Ne_unit;
    //    *Ne = 0;
    //}

    return 1;
}

/*
   Current metric: miknowski spherical coordinates
 */

void gcov_func(double X[NDIM], double gcov[NDIM][NDIM]) {
    MUNULOOP gcov[mu][nu] = 0.;
#if CAR
    gcov[0][0] = -1.;
    gcov[1][1] = 1.;
	gcov[2][2] = 1.;
	gcov[3][3] = 1.;
#endif
    //gcov[2][2] = pow(X[1], 2);
    //gcov[3][3] = pow(X[1] * sin(X[2]), 2);
}

void gcon_func(double gcov[][NDIM], double gcon[][NDIM]) {
    MUNULOOP gcon[mu][nu] = 0.;
#if CAR
	gcon[0][0]=-1.;
	gcon[1][1]=1.;
	gcon[2][2]=1.;
	gcon[3][3]=1.;
#else
    gcon[0][0] = 1 / gcov[0][0];
    gcon[1][1] = 1 / gcov[1][1];
    gcon[2][2] = 1 / (gcov[2][2] + 1e-10);
    gcon[3][3] = 1 / (gcov[3][3] + 1e-10);
#endif
}

/*

   connection calculated analytically for modified Kerr-Schild
    coordinates


   this gives the connection coefficient
   \Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to, e.g., {t,ln(r),theta,phi}
 */

#define DEL (1.e-7)
void get_connection(double X[NDIM], double conn[NDIM][NDIM][NDIM]) {
#if (CAR)
	for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            for (int k = 0; k < NDIM; k++) {
                conn[i][j][k] = 0.;
			}
		}
	}
#else
    double tmp[NDIM][NDIM][NDIM];
    double Xh[NDIM], Xl[NDIM];
    double gcon[NDIM][NDIM];
    double gcov[NDIM][NDIM];
    double gh[NDIM][NDIM];
    double gl[NDIM][NDIM];

    gcov_func(X, gcov);
    gcon_func(gcov, gcon);

    for (int k = 0; k < NDIM; k++) {
        for (int l = 0; l < NDIM; l++)
            Xh[l] = X[l];
        for (int l = 0; l < NDIM; l++)
            Xl[l] = X[l];
        Xh[k] += DEL;
        Xl[k] -= DEL;
        gcov_func(Xh, gh);
        gcov_func(Xl, gl);

        for (int i = 0; i < NDIM; i++) {
            for (int j = 0; j < NDIM; j++) {
                conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
            }
        }
    }

    // Rearrange to find \Gamma_{ijk}
    for (int i = 0; i < NDIM; i++)
        for (int j = 0; j < NDIM; j++)
            for (int k = 0; k < NDIM; k++)
                tmp[i][j][k] =
                    0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);

    // G_{ijk} -> G^i_{jk}
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            for (int k = 0; k < NDIM; k++) {
                conn[i][j][k] = 0.;
                for (int l = 0; l < NDIM; l++)
                    conn[i][j][k] += gcon[i][l] * tmp[l][j][k];
            }
        }
    }
#endif
}
#undef DEL

/* stopping criterion for geodesic integrator */
/* K not referenced intentionally */


#define ROULETTE 1.e3
int stop_criterion(struct of_photon *ph) {
    double wmin = weight_min; /* stop if weight is below minimum weight */
	double r=sqrt( ph->X[1] * ph->X[1] + ph->X[2] * ph->X[2]+ph->X[3]*ph->X[3]);

    if (r> rscreen) {
        if (ph->w < wmin) {
            if (monty_rand() <= 1. / ROULETTE) {
                ph->w *= ROULETTE;
            } else
                ph->w = 0.;
        }
        return 1;
    }

    if (ph->w < wmin) {
        if (monty_rand() <= 1. / ROULETTE) {
            ph->w *= ROULETTE;
        } else {
            ph->w = 0.;
            return 1;
        }
    }
    return (0);
}

/* criterion for recording photon */

int record_criterion(struct of_photon *ph) {
    double r=sqrt( ph->X[1]*ph->X[1] + ph->X[2]*ph->X[2] + ph->X[3]*ph->X[3]);
    if (r > rscreen) return (1);
    else return (0);
}

/* EPS really ought to be related to the number of
   zones in the simulation. */

double clip(double x){
    if (x<STEPMIN) return STEPMIN;
    else if (x>STEPMAX) return STEPMAX;
    else return x;
}
double stepsize(double X[NDIM], double K[NDIM], double tau) {
	double r=sqrt(X[1]*X[1]+X[2]*X[2]);

    if (X[1] < startx[1] || X[1] > stopx[1] || X[2] < startx[2] ||
        X[2] > stopx[2] || X[3] < startx[3] || X[3] > stopx[3]){
            return (r*EPS/(fabs(K[1])+fabs(K[2])+fabs(K[3])+SMALL));
        }
    else return clip(r*EPS)/(fabs(K[1])+fabs(K[2])+fabs(K[3])+SMALL);
}

/*
   record contribution of super photon to spectrum.

   This routine should make minimal assumptions about the
   coordinate system.

 */

void record_super_photon(struct of_photon *ph) {
    int iE, ix2, ix3, ii, jj, kk;
    double lE, dx2, dx3;
    double r, theta, phi;
    if (isnan(ph->w) || isnan(ph->E)) {
        fprintf(stderr, "record isnan: %g %g\n", ph->w, ph->E);
        return;
    }
#pragma omp critical(MAXTAU)
    {
        if (ph->tau_scatt > max_tau_scatt)
            max_tau_scatt = ph->tau_scatt;
    }

#if (CAR)
    r = sqrt(ph->X[1]* ph->X[1]+ ph->X[2]* ph->X[2] + ph->X[3]* ph->X[3]);
    theta = (acos(ph->X[3] / r));
    phi = fmod(atan2(ph->X[2], ph->X[1]), 2 * M_PI);
    if (phi < 0)
        phi += 2 * M_PI;
    if (theta < 0 || theta > M_PI)
        fprintf(stderr, "issue with theta %e\n", theta);
#endif
    /***************************************************/
    /*  get theta bin, without folding around equator   */
    dx2 = (M_PI) / ((double)N_THBINS);
    ix2 = (int)(theta / dx2);
    dx3 = (2 * M_PI) / ((double)N_PHIBINS);
    ix3 = (int)(phi / dx3);
    /***************************************************/

    /* check limits */
    if (ix2 < 0 || ix2 >= N_THBINS) {
        fprintf(stderr, "outside theta bin %e %d\n", theta, ix2);
        return;
    }
    if (ix3 < 0 || ix3 >= N_PHIBINS) {
        fprintf(stderr, "outside phi bin %e %d\n", phi, ix3);
        return;
    }

    /* get energy bin */
    lE = log(ph->E);
    iE =
        (int)((lE - lE0) / dlE + 2.5) - 2; /* bin is centered on iE*dlE + lE0 */
    /* check limits */
    if (iE < 0 || iE >= N_EBINS) {
        return;
    }

#pragma omp atomic
    N_superph_recorded++;
#pragma omp atomic
    N_scatt += ph->nscatt;
    /* sum in photon */
    spect[ix2][ix3][iE].dNdlE += ph->w;
    spect[ix2][ix3][iE].dEdlE += ph->w * ph->E;
    spect[ix2][ix3][iE].tau_abs += ph->w * ph->tau_abs;
    spect[ix2][ix3][iE].tau_scatt += ph->w * ph->tau_scatt;
    spect[ix2][ix3][iE].X1iav += ph->w * ph->X1i;
    spect[ix2][ix3][iE].X2isq += ph->w * (ph->X2i); // * ph->X2i);
    spect[ix2][ix3][iE].X3fsq += ph->w * (ph->X3i); // ph->X3i);
    spect[ix2][ix3][iE].ne0 += ph->w * (ph->ne0);
    spect[ix2][ix3][iE].b0 += ph->w * (ph->b0);
    spect[ix2][ix3][iE].thetae0 += ph->w * (ph->thetae0);
    spect[ix2][ix3][iE].nscatt += ph->nscatt;
    spect[ix2][ix3][iE].nph += 1.;
}

void omp_reduce_spect() {
    /* Combine partial spectra from each OpenMP process		*
     * Inefficient, but only called once so doesn't matter	*/

    int i, j, k;

#pragma omp critical(UPDATE_SPECT)
    {
        for (i = 0; i < N_THBINS; i++) {
            for (k = 0; k < N_PHIBINS; k++) {
                for (j = 0; j < N_EBINS; j++) {
                    shared_spect[i][k][j].dNdlE += spect[i][k][j].dNdlE;
                    shared_spect[i][k][j].dEdlE += spect[i][k][j].dEdlE;
                    shared_spect[i][k][j].tau_abs += spect[i][k][j].tau_abs;
                    shared_spect[i][k][j].tau_scatt += spect[i][k][j].tau_scatt;
                    shared_spect[i][k][j].X1iav += spect[i][k][j].X1iav;
                    shared_spect[i][k][j].X2isq += spect[i][k][j].X2isq;
                    shared_spect[i][k][j].X3fsq += spect[i][k][j].X3fsq;
                    shared_spect[i][k][j].ne0 += spect[i][k][j].ne0;
                    shared_spect[i][k][j].b0 += spect[i][k][j].b0;
                    shared_spect[i][k][j].thetae0 += spect[i][k][j].thetae0;
                    shared_spect[i][k][j].nscatt += spect[i][k][j].nscatt;
                    shared_spect[i][k][j].nph += spect[i][k][j].nph;
                }
            }
        }
    }
#pragma omp barrier
#pragma omp master
    {
        for (i = 0; i < N_THBINS; i++) {
            for (k = 0; k < N_PHIBINS; k++) {
                for (j = 0; j < N_EBINS; j++) {
                    spect[i][k][j].dNdlE = shared_spect[i][k][j].dNdlE;
                    spect[i][k][j].dEdlE = shared_spect[i][k][j].dEdlE;
                    spect[i][k][j].tau_abs = shared_spect[i][k][j].tau_abs;
                    spect[i][k][j].tau_scatt = shared_spect[i][k][j].tau_scatt;
                    spect[i][k][j].X1iav = shared_spect[i][k][j].X1iav;
                    spect[i][k][j].X2isq = shared_spect[i][k][j].X2isq;
                    spect[i][k][j].X3fsq = shared_spect[i][k][j].X3fsq;
                    spect[i][k][j].ne0 = shared_spect[i][k][j].ne0;
                    spect[i][k][j].b0 = shared_spect[i][k][j].b0;
                    spect[i][k][j].thetae0 = shared_spect[i][k][j].thetae0;
                    spect[i][k][j].nscatt = shared_spect[i][k][j].nscatt;
                    spect[i][k][j].nph = shared_spect[i][k][j].nph;
                }
            }
        }
    }
}

/*

   output spectrum to file

 */

#define SPECTRUM_FILE_NAME "grmonty_pwl.spec"

void report_spectrum(double N_superph_made) {
    double JY_FACTOR = 1.e23;
    int i, j, k;
    double dx2, dx3, dOmega, nuLnu, tau_scatt, L, Inu,d;
    FILE *fp;
    char filename[100];

    sprintf(filename, "spec_%d.dat", index);
    fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "trouble opening spectrum file\n");
        exit(0);
    }
#if MPI
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif

    /* output */
    max_tau_scatt = 0.;
    L = 0.;
    for (i = 0; i < N_EBINS; i++) {

        /* output log_10(photon energy/(me c^2)) */
        fprintf(fp, "%10.5g ", exp(i * dlE + lE0) * ME * CL * CL / HPL);

        for (k = 0; k < N_PHIBINS; k++) {
            for (j = 0; j < N_THBINS; j++) {
                /* convert accumulated photon number in each bin
                   to \nu L_\nu, in units of Lsun */
                // to fold the SED


                dx2 = (1.) / (1. * N_THBINS);
                dx3 = (2 * M_PI) / (1. * N_PHIBINS);

                /* factor of 2 accounts for folding around equator */
                // to fold the sed


                dOmega =
                    dOmega_func(j * dx2, (j + 1) * dx2, k * dx3, (k + 1) * dx3);
                d = SOURCE_DIST;

                Inu = 1. / (4. * M_PI * d * d) * JY_FACTOR;
                nuLnu = (ME * CL * CL) * (4. * M_PI / dOmega) * (1. / dlE);
#if MPI
                nuLnu *= shared_spect[j][k][i].dEdlE / ((double)world_size);
#endif

#if (OPENMP && (!MPI))
                nuLnu *= shared_spect[j][k][i].dEdlE;
#endif
                
#if MPI
                tau_scatt =
                    shared_spect[j][k][i].tau_scatt / ((double)world_size) /
                    (shared_spect[j][k][i].dNdlE / ((double)world_size) +
                     SMALL);
#endif

#if (OPENMP && (!MPI))
                tau_scatt = shared_spect[j][k][i].tau_scatt /
                            (shared_spect[j][k][i].dNdlE + SMALL);
#endif
                Inu *= nuLnu / (exp(i * dlE + lE0) * ME * CL * CL / HPL);

                fprintf(fp, "%10.5g %10.5g ", nuLnu, Inu);
                if (tau_scatt > max_tau_scatt)
                    max_tau_scatt = tau_scatt;
                L += nuLnu * dOmega * dlE;
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(stderr,"luminosity %g max_tau_scatt %g\n", L, max_tau_scatt);
    fprintf(stderr, "\n");
    fprintf(stderr, "N_superph_made: %g\n", N_superph_made);
    fprintf(stderr, "N_superph_recorded: %g\n", N_superph_recorded);

    fclose(fp);
}

void sum_struct_ts(void *in, void *inout, int *len, MPI_Datatype *type) {
    /* ignore type, just trust that it's our struct type */

    struct of_spectrum *invals = in;
    struct of_spectrum *inoutvals = inout;
    int i, j, k;

    for (int i = 0; i < *len; i++) {
        inoutvals[i].dNdlE += invals[i].dNdlE;
        inoutvals[i].dEdlE += invals[i].dEdlE;
        inoutvals[i].tau_abs += invals[i].tau_abs;
        inoutvals[i].tau_scatt += invals[i].tau_scatt;
        inoutvals[i].X1iav += invals[i].X1iav;
        inoutvals[i].X2isq += invals[i].X2isq;
        inoutvals[i].X3fsq += invals[i].X3fsq;
        inoutvals[i].ne0 += invals[i].ne0;
        inoutvals[i].b0 += invals[i].b0;
        inoutvals[i].thetae0 += invals[i].thetae0;
        inoutvals[i].nscatt += invals[i].nscatt;
        inoutvals[i].nph += invals[i].nph;
    }
    return;
}

#if MPI
void mpi_reduce_spect() {
    /* Combine partial spectra from each OpenMP process		*
     * Inefficient, but only called once so doesn't matter	*/
    void sum_struct_ts(void *in, void *inout, int *len, MPI_Datatype *type);
    MPI_Datatype tstype;
    MPI_Op sumstruct;

    const int count = 13;
    int blocklens[count];
    MPI_Datatype types[count];
    MPI_Aint disps[count];

    for (int i = 0; i < count; i++) {
        types[i] = MPI_DOUBLE;
        blocklens[i] = 1;
    }

    disps[0] = offsetof(type_spectr, dNdlE);
    disps[1] = offsetof(type_spectr, dEdlE);
    disps[2] = offsetof(type_spectr, nph);
    disps[3] = offsetof(type_spectr, nscatt);
    disps[4] = offsetof(type_spectr, X1iav);
    disps[5] = offsetof(type_spectr, X2isq);
    disps[6] = offsetof(type_spectr, X3fsq);
    disps[7] = offsetof(type_spectr, tau_abs);
    disps[8] = offsetof(type_spectr, tau_scatt);
    disps[9] = offsetof(type_spectr, ne0);
    disps[10] = offsetof(type_spectr, thetae0);
    disps[11] = offsetof(type_spectr, b0);
    disps[12] = offsetof(type_spectr, E0);

    MPI_Type_create_struct(count, blocklens, disps, types, &tstype);
    MPI_Type_commit(&tstype);

    MPI_Op_create(sum_struct_ts, 1, &sumstruct);

    //      MPI_Reduce(&(spect[0][0][0]), &(shared_spect[0][0][0]),
    //      N_THBINS*N_PHIBINS*N_EBINS, tstype, MPI_SUM, 0,
    //               MPI_COMM_WORLD);
    for (int i = 0; i < N_THBINS; i++) {
        for (int j = 0; j < N_PHIBINS; j++) {
            for (int k = 0; k < N_EBINS; k++) {
                MPI_Reduce(&(spect[i][j][k]), &(shared_spect[i][j][k]), 1,
                           tstype, sumstruct, 0, MPI_COMM_WORLD);
            }
        }
    }
}
#endif
