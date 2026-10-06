#include "decs.h"
#include "PLUTO_model.h"
/*

   get HARM simulation data from fname

   checks for consistency of coordinates in data file with
   values of coordinate parameters

   Uses standard HARM data file format

   CFG 1 Sept 07

 */

double dx1min=1.e40,dx2min=1.e40,dx3min=1.e40,dx1max=0.,dx2max=0.,dx3max=0.;
double *x1l,*x1r,*x2l,*x2r,*x3l,*x3r;

void init_pluto_data(char *fname) {
    bias_norm=0.;
    double dV,V=0.;
    int world_rank;
#if MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0)
#endif
		{
	#if (!CAR)
		fprintf(stderr, "\nDon't support metric\n");
		exit(1);
	#endif
        fprintf(stderr, "\nReading HEADER...\n");
		}
	int i,j,k,m,n,success;
    double devi;
    char line[256];
	char* signal;
    int token,re;
    fpos_t file_pos;
    FILE *fp;

    //reading grid.out
    fp = fopen("grid.out","r");
    if (fp == NULL){
        perror("grid.out not found !\n");
        exit(1);
    }

    //Move file pointer to the first line of grid.out that does not begin with a "#".
    success = 0;
    while(!success){
        fgetpos(fp, &file_pos);
        signal=fgets(line,256,fp);
        if (line[0] != '#') success = 1;
        }
    fsetpos(fp, &file_pos);

    //read the left and right sides of x1, x2, x3
    re=fscanf(fp,"%d \n",&(N1));
    x1l=(double *)malloc(N1*sizeof(double));
    x1r=(double *)malloc(N1*sizeof(double));
    dx[1]=(double *)malloc(N1*sizeof(double));
    for(i=0;i<N1;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x1l[i], &x1r[i]);
        dx[1][i]=x1r[i]-x1l[i];
        if(dx1min>dx[1][i]) dx1min=dx[1][i];
        if(dx1max<dx[1][i]) dx1max=dx[1][i];
        }

    re=fscanf(fp,"%d \n",&(N2));
    x2l=(double *)malloc(N2*sizeof(double));
    x2r=(double *)malloc(N2*sizeof(double));
    dx[2]=(double *)malloc(N2*sizeof(double));
    for(i=0;i<N2;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x2l[i], &x2r[i]);
        dx[2][i]=x2r[i]-x2l[i];
        if(dx2min>dx[2][i]) dx2min=dx[2][i];
        if(dx2max<dx[2][i]) dx2max=dx[2][i];
        }

#if TRIM
	int zi,zf,N3g;
	re=fscanf(fp,"%d \n",&(N3g));
	zi=0,zf=N3g-1;

	double *x3l_temp,*x3r_temp;
	x3l_temp=(double*)malloc(N3g*sizeof(double));
    x3r_temp=(double*)malloc(N3g*sizeof(double));

	 for(i=0;i<N3g;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x3l_temp[i], &x3r_temp[i]);
		if ( trim_zi>=x3l_temp[i] && trim_zi<=x3r_temp[i] ) zi=i;
		if ( trim_zf>=x3l_temp[i] && trim_zf<=x3r_temp[i] ) zf=i;
        }
	N3=zf-zi;

	x3l=(double *)malloc(N3*sizeof(double));
    x3r=(double *)malloc(N3*sizeof(double));
    dx[3]=(double *)malloc(N3*sizeof(double));
	for(i=0;i<N3;i++){
		x3l[i]=x3l_temp[i+zi];
		x3r[i]=x3r_temp[i+zi];
		dx[3][i]=x3r[i]-x3l[i];
		if(dx3min>dx[3][i]) dx3min=dx[3][i];
        if(dx3max<dx[3][i]) dx3max=dx[3][i];
		}
    
    devi=(x3r[N3-1]+x3l[0])/2.;//perfrom shift to make z=0 at mid-plane
    for(i=0;i<N3;i++){
        x3l[i]=x3l[i]-devi;
        x3r[i]=x3r[i]-devi;
    }
	#if MPI
    if (world_rank == 0)
	#endif
		{
		fprintf(stderr,"\nTrim: %g-%g rg\n",trim_zi,trim_zf);
		fprintf(stderr,"range %d-%d \n",zi,zf);
		}
    free(x3l_temp);
    free(x3r_temp);

#else
	re=fscanf(fp,"%d \n",&(N3));
    x3l=(double *)malloc(N3*sizeof(double));
    x3r=(double *)malloc(N3*sizeof(double));
    dx[3]=(double *)malloc(N3*sizeof(double));
    for(i=0;i<N3;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x3l[i], &x3r[i]);
        dx[3][i]=x3r[i]-x3l[i];
        if(dx3min>dx[3][i]) dx3min=dx[3][i];
        if(dx3max<dx[3][i]) dx3max=dx[3][i];
        }
    devi=(x3r[N3-1]+x3l[0])/2.;//perfrom shift to make z=0 at mid-plane
    for(i=0;i<N3;i++){
        x3l[i]=x3l[i]-devi;
        x3r[i]=x3r[i]-devi;
    }
#endif
    fclose(fp);

    gam = 4./3.;
    startx[0] = 0.;
    startx[1]=x1l[0];
    startx[2]=x2l[0];
    startx[3]=x3l[0];

    stopx[0] = 1.;
    stopx[1] = x1r[N1-1];
    stopx[2] = x2r[N2-1];
    stopx[3] = x3r[N3-1];

    //read .dbl file
    init_storage();
	fp=fopen(fname, "rb");

#if TRIM
    double *buffer = malloc(N1*N2*N3g*sizeof(double));
#else
    double *buffer = malloc(N1*N2*N3*sizeof(double));
#endif

    for(m=0;m<8;m++){
        if(m==0)n=0;//perform order exchange
        else if(m==7)n=1;
        else n=m+1;
    #if TRIM
        size_t elements_read = fread(buffer, sizeof(double), N3g*N2*N1, fp);

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < N1; i++) {
            for (int j = 0; j < N2; j++) {
                double ***p_n = p[n];
                for (int k = 0; k < N3; k++) {
                    size_t src_index = (k+zi)*N2*N1 + j*N1 + i;
                    p_n[i][j][k] = buffer[src_index];
                    }
                }
            }
   #else
        size_t elements_read = fread(buffer, sizeof(double), N3*N2*N1, fp);

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < N1; i++) {
            for (int j = 0; j < N2; j++) {
                double ***p_n = p[n];
                for (int k = 0; k < N3; k++) {
                    size_t src_index = k*N2*N1 + j*N1 + i;
                    p_n[i][j][k] = buffer[src_index];
                    }
                }
            }
     #endif // TRIM
        }
    free(buffer);
    fclose(fp);

	Thetae_unit = (gam - 1.) * ((MP / ME));
    //calculate Internal energy & bias_norm
    double max_ne_doppler = 0.;
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
                p[UU][i][j][k]=p[UU][i][j][k]/(gam-1.);
                dV=dx[1][i]*dx[2][j]*dx[3][k];

                bias_norm+=dV*pow(p[UU][i][j][k]/p[KRHO][i][j][k]*Thetae_unit,2.);
                V+=dV;
                double v2 = p[U1][i][j][k]*p[U1][i][j][k] +
                            p[U2][i][j][k]*p[U2][i][j][k] +
                            p[U3][i][j][k]*p[U3][i][j][k];
                double bound = p[KRHO][i][j][k]*Ne_unit *
                               (1.+sqrt(v2))/sqrt(1.-v2);
                max_ne_doppler = fmax(max_ne_doppler, bound);
            }
        }
    }
    bias_norm /= V;
    /* Fixed Thomson-depth upper bound across the box, including the largest
       possible fluid-frame path boost. Used only for sampling bias, not opacity. */
    double diagonal = sqrt(pow(stopx[1]-startx[1],2) +
                           pow(stopx[2]-startx[2],2) + pow(stopx[3]-startx[3],2));
    max_tau_scatt = fmax(SMALL, max_ne_doppler * SIGMA_THOMSON * L_unit * diagonal);
    fprintf(stderr,"fixed_sampling_tau_bound: %.12g\n",max_tau_scatt);

//    printf("%d %d %d %lf %lf %lf %lf %lf %lf %lf\n",N1,N2,N3,gam,
//           startx[1],startx[2],startx[3],
//           stopx[1],stopx[2],stopx[3]);
//     exit(1);
//
//    i=250;j=200;k=150;
//    printf("\n %e %e %e %e %e %e %e %e\n", p[KRHO][i][j][k],
//                       p[UU][i][j][k], p[U1][i][j][k], p[U2][i][j][k],
//                       p[U3][i][j][k], p[B1][i][j][k], p[B2][i][j][k],
//                       p[B3][i][j][k]);
//    exit(1);
}
