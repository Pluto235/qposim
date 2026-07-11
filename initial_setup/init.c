/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Relativistic magnetized jet with helical magnetic field periodic boundary conditions
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <math.h>

static void   GetJetValues (double x1, double x2, double x3, double *vj);
const double p0=0.018,a=0.5,Rm=0.25,k=sqrt(10);

#if EOS == TAUB
 static double g_gamma;
#endif


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double vjet[256], vamb[256];
  g_gamma = 4./3.;

  GetJetValues (x1, x2, x3, vjet);



  /* -- set jet values -- */
  if ((x1*x1+x2*x2)<=1.0)
    {
     GetJetValues (x1, x2, x3,v);
    }
  else
    {
    /* -- set ambient values -- */
    v[RHO]=g_inputParam[RHO_A];

    v[VX1] = 0.0;
    v[VX3] = 0.0;
    v[VX2] = 0.0;

    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = 0.0;

    v[PRS] = p0;

    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = 0.0;

    #if NTRACER > 0
    v[TRC] = 0.0;
    #endif
    #ifdef PSI_GLM
    v[PSI_GLM] = 0.0;
    #endif
    }
  //g_smallDensity=1.e-3;
}
/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*
 *
 *
 *********************************************************************** */
{
}


/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3;
  double ***bxs, ***bys, ***bzs;
  double prof, vjet[256], vout[NVAR];

  #ifdef STAGGERED_MHD
  DIM_EXPAND(bxs = d->Vs[BX1s];  ,
           bys = d->Vs[BX2s];  ,
           bzs = d->Vs[BX3s];)
  #endif

  x1  = grid->xgc[IDIR];
  x2  = grid->xgc[JDIR];
  x3  = grid->xgc[KDIR];

#if (GEOMETRY == CYLINDRICAL) || (GEOMETRY == CARTESIAN)
  if (side == X3_BEG)
    {  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER)
      {    /* -- cell-centered boundary conditions -- */
      BOX_LOOP(box,k,j,i)
        {
        GetJetValues (x1[i], x2[j], x3[k], vjet); /* Jet Values */
        NVAR_LOOP(nv) vout[nv] = d->Vc[nv][2*KBEG - k - 1][j][i]; /* Ambient */

        vout[VX3] *= -1.0;
        vout[BX1] *= -1.0;
        vout[BX2] *= -1.0;
        #ifdef PSI_GLM
         vjet[PSI_GLM]  = d->Vc[PSI_GLM][k][JBEG][i] - d->Vc[BX2][k][JBEG][i];
         vout[PSI_GLM] *= -1.0;
        #endif

        prof = (sqrt(x1[i]*x1[i]+x2[j]*x2[j]) <= 1.0);
        NVAR_LOOP(nv) d->Vc[nv][k][j][i] = vout[nv] - (vout[nv] - vjet[nv])*prof;
        }
      }
    }
#endif
}

/* **************************************************************** */
void GetJetValues (double x1, double x2, double x3, double *vj)
/*
 *
 *
 *
 *
 ****************************************************************** */
{
  double lor=g_inputParam[LORENTZ],B0,Bp,r;

  r=sqrt(x1*x1+x2*x2);
  if (r < 1.e-9) r = 1.e-9;

  vj[RHO]=g_inputParam[RHO_J];

  vj[VX1] = 0.0;
  vj[VX2] = 0.0;
  vj[VX3] = sqrt(1.0-1.0/lor/lor);


  vj[PRS]=p0*1.5;

  if (g_inputParam[fB]==0)
    {
    B0=0.1;
    vj[BX1]=0.0;
    vj[BX2]=0.0;
    vj[BX3]=B0;
    }
  else if (g_inputParam[fB]==1)
    {
    B0=0.2;
    vj[BX1]=0.0;
    vj[BX2]=B0*(x1<Rm?(r/Rm):(Rm/r));
    vj[BX3]=0.0;
    }
  else
    {
    B0=0.2;
    vj[BX3]=B0/(1.0+pow(r/a,2.0));
    Bp=k*vj[BX3]*r/a;
    vj[BX1]=Bp*x2/r;
    vj[BX2]=-Bp*x1/r;
    }

  vj[AX1] = 0.0;
  vj[AX2] = 0.0;
  vj[AX3] = 0.0;/*-x1*bpol*/


/* -- set tracer value -- */

  #if NTRACER > 0
  vj[TRC] = 1.0;
  #endif
  #ifdef PSI_GLM
  vj[PSI_GLM] = 0.0;
  #endif
}
