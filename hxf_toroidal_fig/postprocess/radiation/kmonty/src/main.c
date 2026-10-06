
/*

   kmonty, based on grmonty

   Using monte carlo method, estimate spectrum of an appropriately
   scaled axisymmetric GRMHD simulation as a function of
   latitudinal viewing angle.

   Input simulation data is assumed to be in dump format provided by
   HARM code.  Location of input file is, at present, hard coded
   (see init_sim_data.c).

   Nph super-photons are generated in total and then allowed
   to propagate.  They are weighted according to the emissivity.
   The photons are pushed by the geodesic equation.
   Their weight decays according to the local absorption coefficient.
   The photons also scatter with probability related to the local
   scattering opacity.

   The electrons are assumed to have a thermal distribution
   function, and to be at the same temperature as the protons.

   CFG 8-17-06

   Implemented synchrotron sampling, 22 Jan 07

   fixed bugs in tetrad/compton scattering routines, 31 Jan 07

   Implemented normalization for output, 6 Feb 07

   Separated out different synchrotron sampling routines
   into separate files, 8 Mar 07

   fixed bug in energy recording; bug used Kcon[0] rather than
   Kcov[0] as energy, 18 Mar 07

   major reorganization to encapsulate problem-dependent parts 5-6 Nov 07

   added accelerated paticles J. Davelaar 2017-2018

 */

#include "decs.h"
#include "external_dust.h"
#include <time.h>

/* defining declarations for global variables */
struct of_geom **geom;
struct of_spectrum ***shared_spect;
double ***shared_Xi_spec;
double ***shared_ispec;

double F_th[N_ESAMP + 1], F_nth[N_ESAMP + 1], wgt[N_ESAMP + 1];
long int Ns, N_scatt=0, N_superph_recorded=0, N_superph_recorded_total=0;

/* some coordinate parameters */
double a;
double R0, Rin, Rh, Rout, Rms;
double hslope;
double startx[NDIM], stopx[NDIM], dx[NDIM];

double dlE, lE0;
double gam;
double dMsim;
double M_unit, L_unit, T_unit;
double RHO_unit, U_unit, B_unit, Ne_unit, Thetae_unit;
double max_tau_scatt, Ladv, dMact, bias_norm;

double a;
double yhigh;

double *neqpar;
double xprobmin[3], xprobmax[3];
int ng[3], nxlone[3], nleafs;
double *eqpar;

int count_leaf, count_node;
int *forest; //= NULL;
struct block *block_info;
int block_size;
int forest_size;

int n1, n2, n3;
int index;

double Ladv, dMact;

double **Xcoord;
double ***Xgrid;
double ***Xbar;

gsl_rng *r;
gsl_integration_workspace *w;

double MBH,R_LOW,R_HIGH;
double EPS,STEPMIN,STEPMAX,RMAX,SOURCE_DIST;
double trim_zi,trim_zf;
void read_model(char *argv[]);


#pragma omp threadprivate(r)

int main(int argc, char *argv[]) {
    long int N_superph_made_local, N_superph_made;
    int quit_flag;
	int myid;
	int world_size=1,world_rank=0;
    struct of_photon ph;
	time_t currenttime, starttime=time(NULL);
    long int seed_base;
    const char *seed_env;
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s MODEL_FILE PLUTO_DATA_FILE OUTPUT_INDEX\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    read_model(argv);
    /* initialize random number generator */
    sscanf(argv[3], "%d", &index);
    seed_env = getenv("KMONTY_SEED");
    seed_base = seed_env != NULL ? strtol(seed_env, NULL, 10)
                                 : (long int)time(NULL) * index;

#if MPI
    #if OPENMP
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    #else
    MPI_Init(NULL, NULL);
    #endif

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message

    if (world_rank == 0) {
        fprintf(stderr, "\nBOOTING UP kmonty\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0) {
        fprintf(stderr, "\n kmonty running\n\n");
    }

    myid = world_rank;
    init_monty_rand(seed_base + 139 * myid);
    Ns /=(world_size);
#endif

#if OPENMP
#pragma omp parallel private(myid)
    {
        myid = omp_get_thread_num();
		/* Arbitrarily picked initial seed */
	#if MPI
		init_monty_rand(seed_base + 139 * myid +
						1337 * (myid + world_rank * omp_get_num_threads()));
	#else
        init_monty_rand(seed_base + 1476 * myid);
	#endif
    }
#endif

    /* spectral bin parameters */
    dlE = 0.25;        /* bin width */
    lE0 = log(NUMIN*HPL/(ME*CL*CL)); /* location of first bin, in electron rest-mass units */

    /* initialize model data, auxiliary variables */
    init_model(argv);
#if MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    validate_electron_model();

    /** main loop **/
    N_superph_made_local = 0;
    N_superph_recorded = 0;
    N_scatt = 0;
    quit_flag = 0;

#if MPI
    if (world_rank == 0)
#endif
        {
        fprintf(stderr,"\n Ns=%ld \n",Ns*world_size);
        fprintf(stderr, "\nMAIN LOOP\n");
        fflush(stderr);
        }


#if OPENMP
#pragma omp parallel private(ph)
#endif
    {
    while (1) {
    int local_done;
    //generate super photons
    #if OPENMP
        #pragma omp critical(MAKE_SPHOT)
    #endif
        {
        if (!quit_flag) make_super_photon(&ph, &quit_flag);
        local_done = quit_flag;
        }
        if (local_done) break;

        /* push them around */
        int igrid = -1;
        track_super_photon(&ph, &N_superph_recorded, igrid);

        /* step */
        /* give interim reports on rates */
    #if OPENMP
        #pragma omp critical(OUTPUT)
    #endif // OPENMP
        {
        N_superph_made_local += 1;
        if (((int)(N_superph_made_local))%100000==0 && N_superph_made_local>0){
            currenttime = time(NULL);
        #if MPI
        	if (world_rank == 0) {
		        fprintf(stderr, "t= %g s, %g ph/s\n",
		                (double)(currenttime - starttime),
		                (double)N_superph_made_local*world_size/(currenttime - starttime));
        	}
		#else
            fprintf(stderr, "t= %g s, %g ph/s\n",
                    (double)(currenttime - starttime),
                    (double)N_superph_made_local/(currenttime - starttime));
        #endif
            fflush(stderr);
            }
        }
        }
    }

    //record spectrum
#if OPENMP
    #pragma omp barrier
    #pragma omp parallel
    omp_reduce_spect();
#endif

#if MPI
    long int N_superph_made_global     = 0;
    long int N_superph_recorded_global = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_reduce_spect();

    MPI_Allreduce(&N_superph_made_local,&N_superph_made_global,
                  1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(&N_superph_recorded,&N_superph_recorded_global,
                  1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    N_superph_made      = N_superph_made_global;
    N_superph_recorded = N_superph_recorded_global;
#else
    N_superph_made      = N_superph_made_local;
    N_superph_recorded = N_superph_recorded;
#endif

#if (!MPI || defined(MASTER_ONLY))
    if (!MPI || world_rank == 0)
#else
    if (world_rank == 0)
#endif
    {
	currenttime=time(NULL);
	fprintf(stderr,"Time: %g s\n",(double)(currenttime-starttime));

    report_spectrum((double)N_superph_made);
    }

#if MPI
    MPI_Finalize();
#endif
	return 0;
}

// Read model parameters from model.in
void read_model(char *argv[]) {
	double Ntot;
    char temp[100], temp2[100];
    FILE *input;
    char inputfile[100];
	int re;

	sscanf(argv[1], "%s", inputfile);
	input = fopen(inputfile, "r");
	if (input == NULL) {
	    fprintf(stderr, "Can't read file %s! Aborting", inputfile);
	    exit(1);
	}

    // Model parameters
    re=fscanf(input, "%s %s %lf", temp, temp2, &MBH);
    re=fscanf(input, "%s %s %lf", temp, temp2, &SOURCE_DIST);
    re=fscanf(input, "%s %s %lf", temp, temp2, &M_unit);
    re=fscanf(input, "%s %s %lf", temp, temp2, &R_LOW);
    re=fscanf(input, "%s %s %lf", temp, temp2, &R_HIGH);

	re=fscanf(input, "%s %s %lf", temp, temp2, &Ntot);
	Ns=(long int)Ntot;
	re=fscanf(input, "%s %s %lf", temp, temp2, &RMAX);
	re=fscanf(input, "%s %s %lf", temp, temp2, &EPS);
    re=fscanf(input, "%s %s %lf", temp, temp2, &STEPMAX);
    re=fscanf(input, "%s %s %lf", temp, temp2, &STEPMIN);
#if TRIM
	re=fscanf(input, "%s %s %lf", temp, temp2, &trim_zi);
	re=fscanf(input, "%s %s %lf", temp, temp2, &trim_zf);
#endif
    // fprintf(stderr,"%g %g %g %g %g %ld %g %g %g %g\n",
    //         MBH,SOURCE_DIST,M_unit,R_LOW,R_HIGH,Ns,RMAX,EPS,STEPMIN,STEPMAX);
    // to cgs units
    MBH*=MSUN;
    SOURCE_DIST*=KPC2CM;
    
    fclose(input);
}
