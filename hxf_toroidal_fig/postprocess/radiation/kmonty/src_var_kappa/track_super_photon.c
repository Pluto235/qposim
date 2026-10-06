
/*

   main transport subroutine for tracking, absorbing,
   and scattering superphotons

   assumes superphotons do not step out of simulation then back in

 */
#include "decs.h"

#define MAXNSTEP 2000000

void track_super_photon(struct of_photon *ph, long int *N_superph_recorded,
                        int igrid) {
    int bound_flag;
    double dtau_scatt, dtau_abs, dtau;
    double bi, bf;
    double alpha_scatti, alpha_scattf;
    double alpha_absi, alpha_absf;
    double dl, x1;
    double nu, Thetae, Ne, B, theta, sigma, beta;
    struct of_photon php;
    double dtauK, frac;
    double bias = 0.;
    double Xi[NDIM], Ki[NDIM], dKi[NDIM], E0;
    double Gcov[NDIM][NDIM], Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    int nstep = 0;
    int ACCZONE = 0;
    double dx_local = 1e100;
    // Don't track zero-weight photons
    if (ph->w < 1) {
        return;
    }

    // isnan(ph->K[2]) || isnan(ph->K[3]) || ph->w == 0.) {
    /* quality control */
    /*                                if (isnan(ph->X[0]) ||
       isnan(ph->X[1]) ||
       isnan(ph->X[2]) ||
       isnan(ph->X[3]) ||
       isnan(ph->K[0]) ||
       isnan(ph->K[1]) ||
       isnan(ph->K[2]) || isnan(ph->K[3]) ) {
       fprintf(stderr, "track_super_photon: bad input photon.\n");
       fprintf(stderr,
       "X0,X1,X2,X3,K0,K1,K2,K3,w,nscatt: %g %g %g %g %g %g %g %g %g %d\n",
       ph->X[0], ph->X[1], ph->X[2], ph->X[3], ph->K[0],
       ph->K[1], ph->K[2], ph->K[3], ph->w, ph->nscatt);
       return;
       }
     */
    dtauK = 2. * M_PI * L_unit / (ME * CL * CL / HBAR);

    /* Initialize opacities */
    // gcov_func(ph->X, Gcov);
    get_fluid_params(ph->X, Gcov, &Ne, &Thetae, &B, &sigma, &beta, Ucon, Ucov,
                     Bcon, Bcov, &ACCZONE, &dx_local, &igrid);
    // fprintf(stdout,"X=%g %g %g\n",ph->X[1],ph->X[2],ph->X[3]);
	// fprintf(stdout,"Ne=%g Theta_e=%g B=%g sigma=%g beta=%g\n",
	// 		Ne,Thetae,B,sigma,beta);
    theta = get_bk_angle(ph->X, ph->K, Ucov, Bcov, B);
    nu = get_fluid_nu(ph->X, ph->K, Ucov);
    alpha_scatti = alpha_inv_scatt(nu, Thetae, Ne);
    alpha_absi = alpha_inv_abs(nu, Thetae, Ne, B, theta, ACCZONE, sigma, beta);
    bi = bias_func(Thetae, ph->w);

    // fprintf(stderr,"theta=%g nu=%g alpha_scatti=%g alpha_absi=%g bi=%g\n",
    //        theta,nu,alpha_scatti,alpha_absi,bi);

    /* Initialize dK/dlam */
    init_dKdlam(ph->X, ph->K, ph->dKdlam);

    while (!stop_criterion(ph)) {
        /* Save initial position/wave vector */
        Xi[0] = ph->X[0];
        Xi[1] = ph->X[1];
        Xi[2] = ph->X[2];
        Xi[3] = ph->X[3];
        Ki[0] = ph->K[0];
        Ki[1] = ph->K[1];
        Ki[2] = ph->K[2];
        Ki[3] = ph->K[3];
        dKi[0] = ph->dKdlam[0];
        dKi[1] = ph->dKdlam[1];
        dKi[2] = ph->dKdlam[2];
        dKi[3] = ph->dKdlam[3];
        E0 = ph->E0s;
        if (ph->w > 1e75)
            fprintf(stderr, "FAULTY photon here, dtau is %e total tau %e", dtau,
                    ph->tau_abs);
        /* evaluate stepsize */
        dl = stepsize(ph->X, ph->K, dx_local);
        //		if(ph->tau_abs>0.25) dl/=10.;
        /* step the geodesic */
        // push_photon4(ph->X, ph->K, ph->dKdlam, dl);
        push_photon(ph->X, ph->K, ph->dKdlam, dl, &(ph->E0s), 0);
        if (stop_criterion(ph))
            break;
        
        /* allow photon to interact with matter, */
        // gcov_func(ph->X, Gcov);
        get_fluid_params(ph->X, Gcov, &Ne, &Thetae, &B, &sigma, &beta, Ucon,
                         Ucov, Bcon, Bcov, &ACCZONE, &dx_local, &igrid);

        // fprintf(stderr,"theta=%g nu=%g a_scatti=%g a_absi=%g bi=%g\n",
        // theta,nu,alpha_scatti,alpha_absi,bi);

        
        if (alpha_absi > 0. || alpha_scatti > 0. || Ne > 0.) {
            bound_flag = 0;
            if (Ne == 0.)
                bound_flag = 1;
            if (!bound_flag) {
                theta = get_bk_angle(ph->X, ph->K, Ucov, Bcov, B);
                nu = get_fluid_nu(ph->X, ph->K, Ucov);
                if (isnan(nu)) {
                    fprintf(stderr,
                            "isnan nu: track_super_photon dl,E0 %g %g\n", dl,
                            E0);
                    fprintf(stderr, "Xi, %g %g %g %g\n", Xi[0], Xi[1], Xi[2],
                            Xi[3]);
                    fprintf(stderr, "Ki, %g %g %g %g\n", Ki[0], Ki[1], Ki[2],
                            Ki[3]);
                    fprintf(stderr, "dKi, %g %g %g %g\n", dKi[0], dKi[1],
                            dKi[2], dKi[3]);
                    exit(1);
                }
            }

            /* scattering optical depth along step */
            if (bound_flag || nu < 0.) {
                dtau_scatt = 0.5 * alpha_scatti * dtauK * dl;
                dtau_abs = 0.5 * alpha_absi * dtauK * dl;
                alpha_scatti = alpha_absi = 0.;
                bias = 0.;
                bi = 0.;
            } else {
                alpha_scattf = alpha_inv_scatt(nu, Thetae, Ne);
                dtau_scatt = (alpha_scattf + alpha_scatti) * 0.5 * dtauK * dl;
                alpha_scatti = alpha_scattf;
                /* absorption optical depth along step */
                alpha_absf = alpha_inv_abs(nu, Thetae, Ne, B, theta, ACCZONE, sigma, beta);
                dtau_abs = (alpha_absf + alpha_absi) * 0.5 * dtauK * dl;
                alpha_absi = alpha_absf;
                bf = bias_func(Thetae, ph->w);
                bias = (bf + bi) * 0.5;
                bi = bf;
                if (dtau_abs < 0) {
                    fprintf(stderr,
                            "[Error: dtau<0]: alpha=%g theta_e=%g Ne=%g B=%g bi=%g\n",
                            alpha_absf,Thetae,Ne, B,bi);
                }
            }
            if (ph->w > 1e75)
                fprintf(stderr, "FAULTY photon here, dtau is %e total tau %e\n",
                        dtau, ph->tau_abs);
		#if (COMPTON)
            x1 = -log(monty_rand());
            php.w = ph->w / bias;
			//fprintf(stderr,"w=%g bias=%g dtau_scatt=%g\n",ph->w,bias,dtau_scatt);
            if (bias * dtau_scatt > x1 && php.w > weight_min) {
                if (isnan(php.w) || isinf(php.w)) {
                    fprintf(stderr,
                            "w isnan in track_super_photon: Ne, bias, ph->w, "
                            "php.w  %g, %g, %g, %g\n",
                            Ne, bias, ph->w, php.w);
                }

                frac = x1 / (bias * dtau_scatt);

                /* Apply absorption until scattering event */
                dtau_abs *= frac;
                if ((dtau_abs) > 100.)
                    return; /* This photon has been absorbed before scattering
                             */

                dtau_scatt *= frac;
                dtau = dtau_abs + dtau_scatt;
                if (dtau_abs < 1.e-3)
                    ph->w *=
                        (1. - dtau / 24. *
                                  (24. - dtau * (12. - dtau * (4. - dtau))));
                else if (isnan(dtau) || dtau < 0)
                    printf("problem with dtau %e %e %e %e %e\n", dtau, dtau_abs,
                           dtau_scatt, frac, bias);
                else
                    ph->w *= exp(-dtau);

                /* Interpolate position and wave vector to scattering event */
                // push_photon4(Xi, Ki, dKi, dl);
                push_photon(Xi, Ki, dKi, dl * frac, &E0, 0);
                ph->X[0] = Xi[0];
                ph->X[1] = Xi[1];
                ph->X[2] = Xi[2];
                ph->X[3] = Xi[3];
                ph->K[0] = Ki[0];
                ph->K[1] = Ki[1];
                ph->K[2] = Ki[2];
                ph->K[3] = Ki[3];
                ph->dKdlam[0] = dKi[0];
                ph->dKdlam[1] = dKi[1];
                ph->dKdlam[2] = dKi[2];
                ph->dKdlam[3] = dKi[3];
                ph->E0s = E0;

                /* Get plasma parameters at new position */
                // gcov_func(ph->X, Gcov);
                get_fluid_params(ph->X, Gcov, &Ne, &Thetae, &B, &sigma, &beta,
                                 Ucon, Ucov, Bcon, Bcov, &ACCZONE, &dx_local,
                                 &igrid);

                if (Ne > 0.) {
		        scatter_super_photon(ph, &php, Ne, Thetae, B, Ucon, Bcon,
                                         Gcov, ACCZONE, sigma, beta);
                    if (ph->w < 1.e-100) { /* must have been a problem popping k
                                              back onto light cone */
                        return;
                    }
                    track_super_photon(&php, N_superph_recorded, igrid);
                }

                theta = get_bk_angle(ph->X, ph->K, Ucov, Bcov, B);
                nu = get_fluid_nu(ph->X, ph->K, Ucov);
                if (nu < 0.) {
                    alpha_scatti = alpha_absi = 0.;
                } else {
                    alpha_scatti = alpha_inv_scatt(nu, Thetae, Ne);
                    alpha_absi =
                        alpha_inv_abs(nu, Thetae, Ne, B, theta, ACCZONE, sigma, beta);
                }
                bi = bias_func(Thetae, ph->w);

                ph->tau_abs += dtau_abs;
                ph->tau_scatt += dtau_scatt;

            } else 
		#endif
            {
                if ((dtau_abs) > (100.))
                    return; /* This photon has been absorbed */
                ph->tau_abs += dtau_abs;
                ph->tau_scatt += dtau_scatt;
                dtau = dtau_abs + dtau_scatt;
                if (dtau < 1.e-3)
                    ph->w *=
                        (1. - dtau / 24. *
                                  (24. - dtau * (12. - dtau * (4. - dtau))));
                else
                    ph->w *= exp(-dtau);

                if (ph->w > 1e75)
                    fprintf(stderr,
                            "FAULTY photon here, dtau is %e total tau %e", dtau,
                            ph->tau_abs);
            }
        }
        //        printf("%e %e %e %e %e\n",dtau, Thetae, Ne, B ,alpha_absi);
        nstep++;

        /* signs that something's wrong w/ the integration */
        if (nstep > MAXNSTEP) {
            fprintf(
                stderr,
                "X1,X2,X3,K1,K2,K3,bias,dl,steps: %g %g %g %g %g %g %g %g %d\n",
                ph->X[1], ph->X[2], ph->X[3], ph->K[1], ph->K[2], ph->K[3],
                bias, dl, nstep);
            break;
        }
    }
    /* accumulate result in spectrum on escape */
    if (record_criterion(ph) && nstep < MAXNSTEP) {
        record_super_photon(ph);
	#if OPENMP
		#pragma omp atomic
	#endif
        (*N_superph_recorded)++;
    }

    /* done! */
    return;
}
