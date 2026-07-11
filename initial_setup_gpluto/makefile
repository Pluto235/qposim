# *********************************************************
#                  gPLUTO makefile
# *********************************************************

pluto:

ARCH         = linux.nvc++.acc.defs
SRC          = $(PLUTO_DIR)/Src
INCLUDE_DIRS = -I. -I$(SRC)
VPATH        = ./:$(SRC)/New:$(SRC):$(SRC)/Time_Stepping:$(SRC)/States

include $(PLUTO_DIR)/Config/$(ARCH)

include $(SRC)/Templates/default_objects.mk

OBJ     = $(DEFAULT_OBJ)
HEADERS = $(DEFAULT_HEADERS)

include $(SRC)/Math_Tools/makefile

ifeq ($(strip $(PARALLEL)), TRUE)
 CFLAGS += -I$(SRC)/Parallel -DPARALLEL
 include $(SRC)/Parallel/makefile
 ifeq ($(strip $(USE_ASYNC_IO)), TRUE)
  CFLAGS += -DUSE_ASYNC_IO
 endif
endif

ifeq ($(strip $(USE_HDF5)), TRUE)
 CFLAGS += -DUSE_HDF5
 OBJ    += hdf5_io.o
endif

ifeq ($($strip $(USE_PNG)), TRUE)
 CFLAGS += -DUSE_PNG
endif

-include local_make

# ---------------------------------------------------------
#         Problem-specific objects and physics module
# ---------------------------------------------------------

OBJ += ppm_states.o
OBJ += rk_step.o
OBJ += update_stage.o
include $(SRC)/RMHD/makefile

pluto: $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE_DIRS) $<

clean:
	@rm -f *.o
	@echo make clean: done

$(OBJ): $(HEADERS)
