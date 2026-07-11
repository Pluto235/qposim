#define  PHYSICS                        RMHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            5

// Physics-dependent declarations

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   EIGHT_WAVES

// User-defined parameters (labels)

#define  MACH                           0
#define  LORENTZ                        1
#define  RHO_J                          2
#define  RHO_A                          3
#define  FB_MODE                        4

// [Beg] user-defined constants (do not change this line)

#define  SHOCK_FLATTENING               MULTID
#define  FALLBACK_MODE                  FALLBACK_RECOMPUTE

// [End] user-defined constants (do not change this line)
