# Compiler
CC      = mpicc

# Directories
SRCDIR  = src
OBJDIR  = build

# Flags
CFLAGS  = -Wall -O3 -fopenmp -I. -I$(SRCDIR) -I$(SRCDIR)/BHAC -I$(SRCDIR)/sphere -I$(SRCDIR)/PLUTO
LDFLAGS = -lm -lgsl -lgslcblas -fopenmp

# PDC keeps the RaptorP dependencies in a sibling tools directory.  System
# installations continue to work when this prefix is absent.
GSL_PREFIX ?= ../../../../tools/raptorp-deps/prefix/usr
ifneq ($(wildcard $(GSL_PREFIX)/include/gsl/gsl_math.h),)
CFLAGS  += -I$(GSL_PREFIX)/include
LDFLAGS += -L$(GSL_PREFIX)/lib/x86_64-linux-gnu \
           -Wl,-rpath,$(abspath $(GSL_PREFIX)/lib/x86_64-linux-gnu)
endif

# Executable
EXE = kmonty

# ------------------------------------------------------------
# Source files
# ------------------------------------------------------------

SRCS1 = main.c compton.c init_geometry.c tetrads.c geodesics.c \
radiation.c jnu_mixed.c hotcross.c track_super_photon.c \
kappa_sampler.c kappa_sampler_num.c electron_distribution.c external_dust.c \
scatter_super_photon.c \
sphere/sphere_model.c sphere/init_sphere_data.c sphere/sphere_utils.c

SRCS2 = main.c compton.c init_geometry.c tetrads.c geodesics.c \
radiation.c jnu_mixed.c hotcross.c track_super_photon.c \
kappa_sampler.c kappa_sampler_num.c electron_distribution.c external_dust.c scatter_super_photon.c \
BHAC/BHAC_model.c BHAC/BHAC_utils.c BHAC/init_BHAC_data.c

SRCS3 = main.c compton.c init_geometry.c tetrads.c geodesics.c \
radiation.c jnu_mixed.c hotcross.c track_super_photon.c \
kappa_sampler.c kappa_sampler_num.c electron_distribution.c external_dust.c scatter_super_photon.c \
PLUTO/PLUTO_model.c PLUTO/PLUTO_utils.c PLUTO/init_PLUTO_data.c

# Add src/ prefix
SRCS1 := $(addprefix $(SRCDIR)/,$(SRCS1))
SRCS2 := $(addprefix $(SRCDIR)/,$(SRCS2))
SRCS3 := $(addprefix $(SRCDIR)/,$(SRCS3))

# Object files in build/
OBJS1 := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS1))
OBJS2 := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS2))
OBJS3 := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS3))

# ------------------------------------------------------------
# Targets
# ------------------------------------------------------------

all: sphere

sphere: $(OBJDIR) $(OBJS1)
	$(CC) $(CFLAGS) $(OBJS1) $(LDFLAGS) -o $(EXE)

bhac: $(OBJDIR) $(OBJS2)
	$(CC) $(CFLAGS) $(OBJS2) $(LDFLAGS) -o $(EXE)

pluto: $(OBJDIR) $(OBJS3)
	$(CC) $(CFLAGS) $(OBJS3) $(LDFLAGS) -o $(EXE)

pluto-no-compton:
	$(MAKE) pluto OBJDIR=build/no-compton EXE=kmonty_no_compton \
		CFLAGS="$(CFLAGS) -DCOMPTON=0"

# ------------------------------------------------------------
# Compile rule
# ------------------------------------------------------------

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/jnu_mixed.o: $(SRCDIR)/kappa_synch_numeric.h


# ------------------------------------------------------------
# Utility
# ------------------------------------------------------------

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(EXE) kmonty_no_compton
