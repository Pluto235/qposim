////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Three-dimensional relativistic magnetized jet with a helical field.
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

namespace {

constexpr double kAmbientPressure = 0.018;
constexpr double kMagneticRadius  = 0.5;
constexpr double kToroidalScale   = 0.25;
constexpr double kPitch           = 3.16227766016837933199889354443;

#pragma acc routine seq
void GetJetValues(double x1, double x2, double x3, double *vj,
                  const RunConfig &run_config);

} // namespace

////////////////////////////////////////////////////////////////////////
void Init(double *v, double x1, double x2, double x3,
          RunConfig &run_config)
////////////////////////////////////////////////////////////////////////
{
  run_config.gamma = 4.0/3.0;

  const double radius2 = x1*x1 + x2*x2;
  if (radius2 <= 1.0) {
    GetJetValues(x1, x2, x3, v, run_config);
  } else {
    v[RHO] = run_config.inputParam[RHO_A];
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = 0.0;
    v[PRS] = kAmbientPressure;
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
}

////////////////////////////////////////////////////////////////////////
void InitDomain(Data *data, Grid *grid)
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void Analysis(Data *data, Grid *grid)
////////////////////////////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary(const Data *data, RBox *box, int side, Grid *grid)
////////////////////////////////////////////////////////////////////////
{
  int i, j, k, nv;

  if (side != X3_BEG) return;

  const int kbeg = grid->lbeg[KDIR];
  double *x1 = grid->xC[IDIR];
  double *x2 = grid->xC[JDIR];
  double *x3 = grid->xC[KDIR];
  RunConfig &run_config = *(data->run_config);

  #pragma acc parallel loop collapse(3) present(data, box)
  BOX_LOOP(*box, i, j, k) {
    double vjet[NVAR];
    double vout[NVAR];

    GetJetValues(x1[i], x2[j], x3[k], vjet, run_config);

    #pragma acc loop seq
    NVAR_LOOP(nv) {
      vout[nv] = data->Vc(i, j, 2*kbeg - k - 1, nv);
    }

    // Reflect the lower boundary outside the nozzle.
    vout[VX3] *= -1.0;
    vout[BX1] *= -1.0;
    vout[BX2] *= -1.0;
    #ifdef PSI_GLM
    vout[PSI_GLM] *= -1.0;
    #endif

    const double radius2 = x1[i]*x1[i] + x2[j]*x2[j];
    const double nozzle = (radius2 <= 1.0 ? 1.0 : 0.0);

    #pragma acc loop seq
    NVAR_LOOP(nv) {
      data->Vc(i, j, k, nv) = vout[nv] - (vout[nv] - vjet[nv])*nozzle;
    }
  }
}

namespace {

////////////////////////////////////////////////////////////////////////
#pragma acc routine seq
void GetJetValues(double x1, double x2, double x3, double *vj,
                  const RunConfig &run_config)
////////////////////////////////////////////////////////////////////////
{
  const double lorentz = run_config.inputParam[LORENTZ];
  double radius = sqrt(x1*x1 + x2*x2);
  if (radius < 1.e-9) radius = 1.e-9;

  vj[RHO] = run_config.inputParam[RHO_J];
  vj[VX1] = 0.0;
  vj[VX2] = 0.0;
  vj[VX3] = sqrt(1.0 - 1.0/(lorentz*lorentz));
  vj[PRS] = 1.5*kAmbientPressure;

  const int field_mode = static_cast<int>(run_config.inputParam[FB_MODE]);
  if (field_mode == 0) {
    const double b0 = 0.1;
    vj[BX1] = 0.0;
    vj[BX2] = 0.0;
    vj[BX3] = b0;
  } else if (field_mode == 1) {
    const double b0 = 0.2;
    vj[BX1] = 0.0;
    vj[BX2] = b0*(radius < kToroidalScale
                  ? radius/kToroidalScale
                  : kToroidalScale/radius);
    vj[BX3] = 0.0;
  } else {
    const double b0 = 0.2;
    vj[BX3] = b0/(1.0 + pow(radius/kMagneticRadius, 2.0));
    const double bphi = kPitch*vj[BX3]*radius/kMagneticRadius;
    vj[BX1] = bphi*x2/radius;
    vj[BX2] = -bphi*x1/radius;
  }

  vj[AX1] = 0.0;
  vj[AX2] = 0.0;
  vj[AX3] = 0.0;

  #if NTRACER > 0
  vj[TRC] = 1.0;
  #endif
  #ifdef PSI_GLM
  vj[PSI_GLM] = 0.0;
  #endif
}

} // namespace
